import sys
import os
from flask import Flask, render_template, request, jsonify
from pathlib import Path
import threading
from Bio.Seq import Seq
import io
from contextlib import redirect_stdout
import pandas as pd
from werkzeug.utils import secure_filename
import tempfile
from queue import Queue
from tkinter import Tk, filedialog
import time

app = Flask(__name__)
app.config['MAX_CONTENT_LENGTH'] = 100 * 1024 * 1024

current_state = {
    'processing': False,
    'gene_list_file': None,
    'output_file': None,
    'last_message': None,
}

# Thread-safe queues for dialog communication
dialog_request_queue = Queue()
dialog_response_queue = Queue()
dialog_worker_ready = threading.Event()

# Import functions from C_Create_Codebook
sys.path.insert(0, os.path.dirname(__file__))
from C_Create_Codebook import (
    get_mouse_sequences_file,
    load_primers_from_file,
    load_readout_sequences_from_file,
    load_raefish_adapter_sequences,
    load_raefish_readout_probes,
    load_raefish_sequences,
    load_gene_list,
    generate_hybridization_codebook
)

def tkinter_dialog_worker():
    """Worker thread that handles all tkinter dialogs"""
    try:
        root = Tk()
        root.withdraw()  # Hide the root window
        root.attributes('-topmost', True)  # Keep window on top
        print("✓ Tkinter dialog worker initialized")
        
        # Signal that the worker is ready
        dialog_worker_ready.set()
        
        while True:
            try:
                # Wait for a dialog request with a timeout
                request_data = dialog_request_queue.get(timeout=1)
                
                if request_data is None:  # Signal to stop
                    break
                
                dialog_type = request_data.get('type')
                result = None
                error = None
                
                if dialog_type == 'select_folder':
                    try:
                        print("Opening folder selection dialog...")
                        result = filedialog.askdirectory(
                            title="Select Output Folder"
                        )
                        if result:
                            print(f"✓ Folder selected: {result}")
                        else:
                            print("User cancelled folder selection")
                            error = "No folder selected"
                    except Exception as e:
                        print(f"✗ Error in folder dialog: {str(e)}")
                        error = str(e)
                
                # Send response back
                dialog_response_queue.put({'result': result, 'error': error})
                
            except Exception as e:
                # Check if this is just a timeout (empty exception from queue.get)
                error_str = str(e)
                if error_str and error_str != "":
                    print(f"✗ Queue error: {error_str}")
                    try:
                        dialog_response_queue.put({'result': None, 'error': error_str})
                    except:
                        pass
                # If it's just a timeout, continue waiting (don't print anything)
            
            # Give other threads a chance to run
            time.sleep(0.01)
        
        root.destroy()
        print("✓ Tkinter dialog worker stopped")
        
    except Exception as e:
        print(f"✗ Fatal error in tkinter worker: {str(e)}")
        dialog_worker_ready.set()

@app.route('/api/upload-gene-list', methods=['POST'])
def upload_gene_list():
    """Handle gene list file upload"""
    try:
        if 'file' not in request.files:
            return jsonify({'success': False, 'error': 'No file provided'}), 400
        
        file = request.files['file']
        if file.filename == '':
            return jsonify({'success': False, 'error': 'No file selected'}), 400
        
        if not file.filename.endswith(('.xlsx', '.xls')):
            return jsonify({'success': False, 'error': 'File must be Excel format (.xlsx or .xls)'}), 400
        
        # Save to temp directory
        temp_dir = tempfile.gettempdir()
        filename = secure_filename(file.filename)
        file_path = os.path.join(temp_dir, filename)
        file.save(file_path)
        
        return jsonify({'success': True, 'file_path': file_path})
    
    except Exception as e:
        return jsonify({'success': False, 'error': str(e)}), 500

@app.route('/api/select-output-folder', methods=['POST'])
def select_output_folder():
    """Open a folder picker dialog using thread-safe queue"""
    try:
        # Check if worker is ready
        if not dialog_worker_ready.is_set():
            print("✗ Dialog worker not ready yet")
            return jsonify({'success': False, 'error': 'Dialog system initializing. Please try again in a moment.'}), 503
        
        print("Sending folder selection request to worker thread...")
        # Send request to tkinter worker thread
        dialog_request_queue.put({'type': 'select_folder'})
        
        # Wait for response with timeout
        try:
            response = dialog_response_queue.get(timeout=120)
        except Exception as e:
            print(f"✗ Timeout or error waiting for dialog: {str(e)}")
            return jsonify({'success': False, 'error': 'Dialog timeout - the folder picker did not respond'}), 504
        
        folder = response.get('result')
        error = response.get('error')
        
        if error:
            print(f"✗ Dialog returned error: {error}")
            return jsonify({'success': False, 'error': error}), 400
        elif folder and folder.strip():
            print(f"✓ Returning folder path to client: {folder}")
            return jsonify({'success': True, 'folder': folder})
        else:
            print("✗ No valid folder path returned")
            return jsonify({'success': False, 'error': 'No folder selected'}), 400
    
    except Exception as e:
        print(f"✗ Exception in select_output_folder: {str(e)}")
        return jsonify({'success': False, 'error': f'Dialog error: {str(e)}'}), 500

@app.route('/')
def index():
    return render_template('index.html')

@app.route('/api/generate-codebook', methods=['POST'])
def generate_codebook():
    if current_state['processing']:
        return jsonify({'success': False, 'error': 'Already processing. Please wait for the current generation to complete.'}), 400

    try:
        data = request.json
        
        # Validate inputs
        mode = data.get('mode')
        gene_count = data.get('gene_count')
        gene_list_file = data.get('gene_list_file')
        hybridizations = data.get('hybridizations', 4)
        channels = data.get('channels', 3)
        output_file = data.get('output_file')
        reverse_complement = data.get('reverse_complement', False)
        barcode_type = data.get('barcode_type', 'combinatorial')
        
        if not mode:
            return jsonify({'success': False, 'error': 'Mode must be selected'}), 400
        
        if mode not in ['MERFISH', 'RAEFISH']:
            return jsonify({'success': False, 'error': 'Invalid mode selected'}), 400
        
        if barcode_type not in ['combinatorial', 'binary']:
            return jsonify({'success': False, 'error': 'Invalid barcode type selected'}), 400
        
        if mode == 'MERFISH':
            # MERFISH can use either gene count or gene list file
            if gene_count is None and (not gene_list_file or not os.path.exists(gene_list_file)):
                return jsonify({'success': False, 'error': 'Either gene count or gene list file must be provided'}), 400
        elif mode == 'RAEFISH':
            # RAEFISH requires gene list file
            if not gene_list_file or not os.path.exists(gene_list_file):
                return jsonify({'success': False, 'error': 'Gene list file not selected or file does not exist. Please select a valid Excel file.'}), 400
        
        if not output_file:
            return jsonify({'success': False, 'error': 'Output file path not specified. Please specify where to save the codebook.'}), 400
        
        # Ensure output directory exists
        output_dir = os.path.dirname(output_file)
        if output_dir and not os.path.exists(output_dir):
            try:
                os.makedirs(output_dir)
            except Exception as e:
                return jsonify({'success': False, 'error': f'Cannot create output directory: {str(e)}'}), 400
        
        try:
            hybridizations = int(hybridizations)
            channels = int(channels)
            if hybridizations < 1 or channels < 1:
                return jsonify({'success': False, 'error': 'Hybridizations and Channels must be positive integers'}), 400
        except ValueError:
            return jsonify({'success': False, 'error': 'Hybridizations and Channels must be valid integers'}), 400
        
        current_state['processing'] = True
        current_state['last_message'] = None
        current_state['output_file'] = output_file
        
        # Run codebook generation in background thread
        def run_generation():
            output_buffer = io.StringIO()
            
            try:
                # Capture stdout
                with redirect_stdout(output_buffer):
                    mouse_sequences_file = get_mouse_sequences_file()
                    READS = 4
                    
                    if mode == "MERFISH":
                        print("Loading MERFISH mode...")
                        readout_sequences_list = load_readout_sequences_from_file(mouse_sequences_file, sheet_name='MERFISH_readout_sequences')
                        five_primer, three_primer = load_primers_from_file(mouse_sequences_file, primer_type='MERFISH')
                        
                        if gene_count is not None:
                            # Use synthetic gene names
                            genes_list = [f"Gene_{i+1}" for i in range(gene_count)]
                            gene_expressions = [1.0] * gene_count
                            print(f"Generating {gene_count} synthetic genes")
                        else:
                            # Load from file
                            genes_list, gene_expressions = load_gene_list(gene_list_file)
                            print(f"Successfully loaded {len(genes_list)} genes from gene list file")
                        
                        generate_hybridization_codebook(mode, len(genes_list), READS, hybridizations, channels, readout_sequences_list, genes_list, gene_expressions, five_primer, three_primer, output_file, reverse_complement=reverse_complement, optimize_balance=data.get('optimize_balance', True), enforce_pairwise=data.get('enforce_pairwise', True), barcode_type=barcode_type)
                        
                    else:  # RAEFISH
                        print("Loading RAEFISH mode...")
                        genes_list, gene_expressions = load_gene_list(gene_list_file)
                        print(f"Successfully loaded {len(genes_list)} genes from gene list file")
                        five_primer, three_primer = load_primers_from_file(mouse_sequences_file, primer_type='RAEFISH')
                        adapter_seqs = load_raefish_adapter_sequences(mouse_sequences_file)
                        adapter_sequences_list = [str(Seq(seq).reverse_complement()) for seq in adapter_seqs]
                        raefish_probes = load_raefish_readout_probes(mouse_sequences_file)
                        raefish_gene_mapping = load_raefish_sequences(mouse_sequences_file)
                        generate_hybridization_codebook(mode, len(genes_list), READS, hybridizations, channels, adapter_sequences_list, genes_list, gene_expressions, five_primer, three_primer, output_file, adapter_seqs, raefish_probes, raefish_gene_mapping, readout_label="Adapter", optimize_balance=data.get('optimize_balance', True), enforce_pairwise=data.get('enforce_pairwise', True), barcode_type=barcode_type)
                
                output_text = output_buffer.getvalue()
                current_state['last_message'] = output_text
                
                # Check if file was generated
                if os.path.exists(output_file):
                    file_size = os.path.getsize(output_file) / (1024 * 1024)  # Size in MB
                    print(f"✓ Codebook generated successfully!")
                    print(f"✓ File saved to: {output_file}")
                    print(f"✓ File size: {file_size:.2f} MB")
                else:
                    print(f"✗ Output file was not created at {output_file}")
                
                current_state['processing'] = False
                
            except FileNotFoundError as e:
                error_msg = f"File not found error: {str(e)}\n\nPlease ensure all required files exist."
                current_state['last_message'] = error_msg
                current_state['processing'] = False
                print(f"✗ Error: {error_msg}")
            
            except PermissionError as e:
                error_msg = f"Permission denied: Cannot save the codebook.\n\n✗ The file may be open in another program (close it first)\n✗ You may not have write permissions in the chosen folder\n✗ The chosen folder may be read-only\n\nPlease check these and try again with a different folder or file name."
                current_state['last_message'] = error_msg
                current_state['processing'] = False
                print(f"✗ Error: {error_msg}")
                
            except ValueError as e:
                error_msg = f"Input validation error: {str(e)}\n\nPlease check your input parameters:\n- Number of genes must be compatible with hybridizations × channels\n- Gene list file must have correct format"
                current_state['last_message'] = error_msg
                current_state['processing'] = False
                print(f"✗ Error: {error_msg}")
                
            except KeyError as e:
                error_msg = f"Column or data error: {str(e)}\n\nThe input file or MOUSE_sequences.xlsx may have an incorrect format."
                current_state['last_message'] = error_msg
                current_state['processing'] = False
                print(f"✗ Error: {error_msg}")
                
            except Exception as e:
                error_msg = f"An unexpected error occurred: {str(e)}\n\nPlease check the console for more details."
                current_state['last_message'] = error_msg
                current_state['processing'] = False
                print(f"✗ Error: {error_msg}")
                import traceback
                traceback.print_exc()
        
        thread = threading.Thread(target=run_generation)
        thread.daemon = True
        thread.start()
        
        return jsonify({'success': True, 'message': 'Codebook generation started in background'})
    
    except Exception as e:
        current_state['processing'] = False
        error_msg = f"Unexpected error: {str(e)}"
        print(f"✗ {error_msg}")
        import traceback
        traceback.print_exc()
        return jsonify({'success': False, 'error': error_msg}), 500

@app.route('/api/check-status', methods=['GET'])
def check_status():
    """Check if generation is complete and return message"""
    return jsonify({
        'processing': current_state['processing'],
        'message': current_state['last_message'],
        'output_file': current_state['output_file']
    })

if __name__ == '__main__':
    print("Starting seqFISH Codebook Generator...")
    print("Initializing dialog system...")
    
    # Start the tkinter dialog worker thread BEFORE Flask
    dialog_worker_thread = threading.Thread(target=tkinter_dialog_worker, daemon=True)
    dialog_worker_thread.start()
    
    # Wait for the worker to be ready (with timeout)
    if dialog_worker_ready.wait(timeout=5):
        print("✓ Dialog system initialized successfully")
    else:
        print("⚠ Warning: Dialog system may not be ready yet")
    
    print("Starting Flask web server...")
    print("Opening browser to http://localhost:5000")
    app.run(debug=False, host='127.0.0.1', port=5000, use_reloader=False, threaded=True)