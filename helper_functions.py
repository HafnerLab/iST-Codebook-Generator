import pandas as pd
import numpy as np
import random
from itertools import combinations
from Bio.Seq import Seq
from pathlib import Path

# CONFIGURATION (only used when running as main script):
MODE = "RAEFISH"
GENES = 16
HYBRIDIZATIONS = 4
CHANNELS = 3
GENE_LIST_FILE = "Y:\\seqFISH+\\Testing_design\\Gene_list_pilot.xlsx"
OUTPUT_FILE = "Y:\\seqFISH+\\Testing_design\\Codebook_output_test_raefish.xlsx"
READS = 4

# Automatically locate MOUSE_SEQUENCES_FILE from doc folder
def get_mouse_sequences_file():
    """Automatically locate MOUSE_sequences.xlsx from the doc folder in the same repository"""
    script_dir = Path(__file__).parent
    doc_folder = script_dir / "doc"
    mouse_sequences_file = doc_folder / "MOUSE_sequences.xlsx"
    
    if not mouse_sequences_file.exists():
        raise FileNotFoundError(f"MOUSE_sequences.xlsx not found in {doc_folder}. "
                               f"Please ensure the file is located at: {mouse_sequences_file}")
    
    return str(mouse_sequences_file)

# Load primer sequences automatically from MOUSE_sequences file
def load_primers_from_file(file_path, primer_type='MERFISH'):
    """Load FIVE_PRIMER and THREE_PRIMER from MOUSE_sequences file"""
    if primer_type == 'MERFISH':
        df_primers = pd.read_excel(file_path, sheet_name='MERFISH_primers', header=0)
        FIVE_PRIMER = df_primers.iloc[0, 2]
        THREE_PRIMER = df_primers.iloc[1, 2]
    elif primer_type == 'RAEFISH':
        df_primers = pd.read_excel(file_path, sheet_name='RAEFISH_primers', header=0)
        THREE_PRIMER = df_primers.iloc[8, 3]
        FIVE_PRIMER = df_primers.iloc[9, 3]
    
    return FIVE_PRIMER, THREE_PRIMER

# Load readout sequences automatically from MOUSE_sequences file
def load_readout_sequences_from_file(file_path, sheet_name='MERFISH_readout_sequences'):
    """Load readout sequences from MOUSE_sequences file"""
    df_readout = pd.read_excel(file_path, sheet_name=sheet_name, header=0)
    readout_sequences = df_readout.iloc[:, 1].dropna().tolist()
    return readout_sequences

# Load RAEFISH adapter sequences
def load_raefish_adapter_sequences(file_path):
    """Load RAEFISH adapter sequences from MOUSE_sequences file"""
    df_adapters = pd.read_excel(file_path, sheet_name='RAEFISH_adapter_sequences', header=0)
    adapter_sequences = df_adapters.iloc[:, 1].dropna().tolist()
    return adapter_sequences

# Load RAEFISH readout probes
def load_raefish_readout_probes(file_path):
    """Load RAEFISH readout probes with landing sequences from MOUSE_sequences file"""
    df_readout_probes = pd.read_excel(file_path, sheet_name='RAEFISH_readout_probes', header=0)
    readout_probes = {}
    for idx, row in df_readout_probes.iterrows():
        channel = int(row.iloc[0])
        landing_sequence = row.iloc[1]
        readout_probes[channel] = landing_sequence
    return readout_probes

# Load RAEFISH sequences and gene mapping
def load_raefish_sequences(file_path):
    """Load RAEFISH sequences and gene mapping from MOUSE_sequences file"""
    df_raefish = pd.read_excel(file_path, sheet_name='RAEFISH_sequences', header=0)
    # Assuming: Column 0 = Ensembl_GeneID, Column 1 = Gene_ShortName, Column 2 = Target_1, Column 3 = Target_2
    # Column 4 = Padlock_1, Column 5 = Padlock_2, Column 6 = Priming_1, Column 7 = Priming_2
    gene_mapping = {}
    for idx, row in df_raefish.iterrows():
        gene_short_name = row.iloc[1]
        ensembl_gene_id = row.iloc[0]
        target_1 = row.iloc[2]
        target_2 = row.iloc[3]
        padlock_1 = row.iloc[4]
        padlock_2 = row.iloc[5]
        priming_1 = row.iloc[6]
        priming_2 = row.iloc[7]
        gene_mapping[gene_short_name] = {
            'Target_1': target_1,
            'Target_2': target_2,
            'Padlock_1': padlock_1,
            'Padlock_2': padlock_2,
            'Priming_1': priming_1,
            'Priming_2': priming_2
        }
        # Also map by Ensembl_GeneID for alternative lookup
        gene_mapping[ensembl_gene_id] = {
            'Target_1': target_1,
            'Target_2': target_2,
            'Padlock_1': padlock_1,
            'Padlock_2': padlock_2,
            'Priming_1': priming_1,
            'Priming_2': priming_2
        }
    return gene_mapping

# Load gene list with expression values
def load_gene_list(file_path):
    """Load gene list with expression values from an Excel file"""
    df_genes = pd.read_excel(file_path, header=None)
    genes = df_genes.iloc[:, 0].tolist()
    
    # Check if expression values are provided (second column)
    expressions = None
    if df_genes.shape[1] > 1:
        try:
            expressions = df_genes.iloc[:, 1].astype(float).tolist()
            # Normalize expression values to 0-1 range
            min_expr = min(expressions)
            max_expr = max(expressions)
            if max_expr > min_expr:
                expressions = [(expr - min_expr) / (max_expr - min_expr) for expr in expressions]
            else:
                expressions = [1.0] * len(expressions)
        except:
            expressions = None
    
    if expressions is None:
        expressions = [1.0] * len(genes)
    
    return genes, expressions


def _build_merfish_to_order_idt(mouse_sequences_file, pattern_df):
    """Build To_order_IDT for MERFISH: primers + channel-prefixed readout probes."""
    channel_prefixes = {1: '/5TYE665N/', 2: '/5Cy3/', 3: '/5Alex488N/'}
    rows = []

    # Primers from MERFISH_primers sheet (col 0 = name, col 2 = sequence)
    df_primers = pd.read_excel(mouse_sequences_file, sheet_name='MERFISH_primers', header=0)
    for _, row in df_primers.iterrows():
        name = str(row.iloc[0]).strip() if pd.notna(row.iloc[0]) else ''
        seq = str(row.iloc[2]).strip() if pd.notna(row.iloc[2]) else ''
        if name and seq and seq.lower() != 'nan':
            rows.append({'Name': name, 'Sequence': seq})

    # Readout probes: numbered, prefixed by channel
    for i, (_, row) in enumerate(pattern_df.iterrows()):
        channel = int(row['Channel'])
        prefix = channel_prefixes.get(channel, '')
        seq = str(row['Readout_sequence']).strip()
        rows.append({'Name': f'Readout_{i + 1}', 'Sequence': f'{prefix}{seq}'})

    return pd.DataFrame(rows, columns=['Name', 'Sequence'])


def _build_raefish_to_order_twist(genes, encoding_probes_rc_list1, encoding_probes_rc_list2,
                                   padlock_1_list, padlock_2_list, priming_1_list, priming_2_list):
    """Build To_order_Twist for RAEFISH: all probes flattened to one row per probe."""
    import math
    probe_entries = [
        ('Encoding_probe1_RC', encoding_probes_rc_list1),
        ('Encoding_probe2_RC', encoding_probes_rc_list2),
        ('Padlock_1', padlock_1_list),
        ('Padlock_2', padlock_2_list),
        ('Priming_1', priming_1_list),
        ('Priming_2', priming_2_list),
    ]
    rows = []
    for probe_type, probe_list in probe_entries:
        for gene, seq in zip(genes, probe_list):
            if seq is not None and not (isinstance(seq, float) and math.isnan(seq)):
                rows.append({'Name': f'{gene}_{probe_type}', 'Sequence': str(seq)})
    return pd.DataFrame(rows, columns=['Name', 'Sequence'])


def _build_raefish_to_order_idt(mouse_sequences_file, pattern_df_raefish):
    """Build To_order_IDT for RAEFISH: primers + adapter probes + readout/blocking oligos."""
    import math
    rows = []

    # Primers from RAEFISH_primers: combine col 0 and col 1 for name, col 2 for sequence
    df_primers = pd.read_excel(mouse_sequences_file, sheet_name='RAEFISH_primers', header=0)
    for _, row in df_primers.iterrows():
        col0 = str(row.iloc[0]).strip() if pd.notna(row.iloc[0]) else ''
        col1 = str(row.iloc[1]).strip() if pd.notna(row.iloc[1]) else ''
        name = f'{col0} {col1}'.strip()
        seq = str(row.iloc[2]).strip() if pd.notna(row.iloc[2]) else ''
        if name and seq and seq.lower() != 'nan':
            rows.append({'Name': name, 'Sequence': seq})

    # Adapter probes from Sequence_Pattern
    if 'Adapter_probe' in pattern_df_raefish.columns:
        for i, (_, row) in enumerate(pattern_df_raefish.iterrows()):
            probe = row['Adapter_probe']
            if probe is not None and not (isinstance(probe, float) and math.isnan(probe)):
                rows.append({'Name': f'Adapter_probe_{i + 1}', 'Sequence': str(probe)})

    # Readout probes and blocking oligos from RAEFISH_readout_probes
    # col 0 (Channel) = channel number, col 1 (Readout_probes) = sequence, col 3 (Blocking_oligos) = sequence
    df_readout = pd.read_excel(mouse_sequences_file, sheet_name='RAEFISH_readout_probes', header=0)
    for _, row in df_readout.iterrows():
        channel = row.iloc[0]
        val = row.iloc[1]
        if val is not None and not (isinstance(val, float) and math.isnan(val)) and str(val).strip():
            rows.append({'Name': f'Readout_Ch{int(channel)}', 'Sequence': str(val).strip()})
        blocking = row.iloc[3]
        if blocking is not None and not (isinstance(blocking, float) and math.isnan(blocking)) and str(blocking).strip():
            rows.append({'Name': f'Blocking_oligo_Ch{int(channel)}', 'Sequence': str(blocking).strip()})

    return pd.DataFrame(rows, columns=['Name', 'Sequence'])


def generate_hybridization_codebook(MODE, GENES, READS, HYBRIDIZATIONS, CHANNELS, READOUT_SEQUENCES_LIST, genes, gene_expressions, FIVE_PRIMER, THREE_PRIMER, OUTPUT_FILE, raefish_adapter_sequences=None, raefish_readout_probes=None, raefish_gene_mapping=None, readout_label="Readout", reverse_complement=False, optimize_balance=True, enforce_pairwise=True, barcode_type="combinatorial", mouse_sequences_file=None):
    """Generate hybridization codebook and save to Excel
    
    Parameters:
    -----------
    reverse_complement : bool, default=False
        If True, applies reverse complement to five_sec_seq and three_sec_seq.
        If False, leaves sequences intact (for PaintSHOP).
        Only applicable for MERFISH mode.
    
    optimize_balance : bool, default=True
        If True, optimizes the frequency/expression balance across readouts using iterative selection.
        If False, simply maximizes the number of valid barcodes without iterative optimization.
        In both cases, the Hamming distance constraint (at least 2 positions different) is maintained.
    
    enforce_pairwise : bool, default=True
        If True, enforces that no two barcodes' union appears in pair_set (strict uniqueness constraint).
        If False, only enforces Hamming distance constraint (at least 2 positions different).
    
    barcode_type : str, default="combinatorial"
        Type of barcode to generate:
        - "combinatorial": Barcodes as combinations of readout numbers (current behavior)
        - "binary": Barcodes as binary bit codes (0s and 1s), one bit per imaging round
    """
    mouse_seq_file = mouse_sequences_file if mouse_sequences_file is not None else get_mouse_sequences_file()
    strings = READOUT_SEQUENCES_LIST
    max_rows = len(strings)

    # Calculate the minimum sequences needed
    min_sec_seq = HYBRIDIZATIONS * CHANNELS
    if min_sec_seq > max_rows:
        raise ValueError(f"The readout sequences must contain at least {min_sec_seq} sequences, but it has only {max_rows} sequences.")

    n_numbers = HYBRIDIZATIONS * CHANNELS

    # Step 2: make codebook
    def create_simple_codebook(READS, n_numbers, GENES):
        """Simple method: maximize valid barcodes without iterative optimization"""
        all_combinations = list(combinations(range(1, n_numbers + 1), READS))
        random.shuffle(all_combinations)

        valid_combinations = []
        pair_set = set() if enforce_pairwise else None

        for combo in all_combinations:
            # Always check Hamming distance constraint: at least 2 positions must differ
            is_valid = all(len(set(combo) ^ set(existing)) >= 2 for existing in valid_combinations)

            # Only check pairwise uniqueness constraint if enforce_pairwise is True
            if is_valid and enforce_pairwise:
                for existing_combo in valid_combinations:
                    combined_set = frozenset(set(combo) | set(existing_combo))
                    if combined_set in pair_set:
                        is_valid = False
                        break

            if is_valid:
                valid_combinations.append(combo)
                # Only update pair_set if enforce_pairwise is True
                if enforce_pairwise:
                    for existing_combo in valid_combinations:
                        pair_set.add(frozenset(set(combo) | set(existing_combo)))

        if len(valid_combinations) < GENES:
            raise ValueError(f"Unable to generate {GENES} unique codebook combinations with the constraints. "
                            f"Generated only {len(valid_combinations)} valid combinations.")

        # Simply take the first GENES valid combinations without optimization
        simple_codebook = valid_combinations[:GENES]
        
        frequencies = np.zeros(n_numbers, dtype=int)
        for combo in simple_codebook:
            for num in combo:
                frequencies[num - 1] += 1

        print(f"Frequencies of each {readout_label} in the codebook:")
        for i in range(1, n_numbers + 1):
            print(f"{readout_label} {i}: {frequencies[i - 1]}")

        codebook_df = pd.DataFrame(simple_codebook, columns=[f'Read{i+1}' for i in range(READS)])
        return codebook_df
    
    def create_max_codebook(READS, n_numbers, GENES):
        """Original method: balance frequencies only"""
        all_combinations = list(combinations(range(1, n_numbers + 1), READS))
        random.shuffle(all_combinations)

        valid_combinations = []
        pair_set = set() if enforce_pairwise else None

        for combo in all_combinations:
            # Always check Hamming distance constraint: at least 2 positions must differ
            is_valid = all(len(set(combo) ^ set(existing)) >= 2 for existing in valid_combinations)

            # Only check pairwise uniqueness constraint if enforce_pairwise is True
            if is_valid and enforce_pairwise:
                for existing_combo in valid_combinations:
                    combined_set = frozenset(set(combo) | set(existing_combo))
                    if combined_set in pair_set:
                        is_valid = False
                        break

            if is_valid:
                valid_combinations.append(combo)
                # Only update pair_set if enforce_pairwise is True
                if enforce_pairwise:
                    for existing_combo in valid_combinations:
                        pair_set.add(frozenset(set(combo) | set(existing_combo)))

        if len(valid_combinations) < GENES:
            raise ValueError(f"Unable to generate {GENES} unique codebook combinations with the constraints. "
                            f"Generated only {len(valid_combinations)} valid combinations.")

        max_codebook = []
        max_balance = float('inf')

        for _ in range(10000):
            random.shuffle(valid_combinations)
            temp_codebook = valid_combinations[:GENES]
            
            temp_frequencies = np.zeros(n_numbers, dtype=int)
            for combo in temp_codebook:
                for num in combo:
                    temp_frequencies[num - 1] += 1
            
            balance_score = np.std(temp_frequencies)
            if balance_score < max_balance:
                max_balance = balance_score
                max_codebook = temp_codebook
                frequencies = temp_frequencies

        print(f"Frequencies of each {readout_label} in the codebook:")
        for i in range(1, n_numbers + 1):
            print(f"{readout_label} {i}: {frequencies[i - 1]}")

        codebook_df = pd.DataFrame(max_codebook, columns=[f'Read{i+1}' for i in range(READS)])
        return codebook_df

    def create_balanced_codebook(READS, n_numbers, GENES, gene_expressions):
        """Expression-based method: balance expression levels across numbers"""
        all_combinations = list(combinations(range(1, n_numbers + 1), READS))
        random.shuffle(all_combinations)

        valid_combinations = []
        pair_set = set() if enforce_pairwise else None

        for combo in all_combinations:
            # Always check Hamming distance constraint: at least 2 positions must differ
            is_valid = all(len(set(combo) ^ set(existing)) >= 2 for existing in valid_combinations)

            # Only check pairwise uniqueness constraint if enforce_pairwise is True
            if is_valid and enforce_pairwise:
                for existing_combo in valid_combinations:
                    combined_set = frozenset(set(combo) | set(existing_combo))
                    if combined_set in pair_set:
                        is_valid = False
                        break

            if is_valid:
                valid_combinations.append(combo)
                # Only update pair_set if enforce_pairwise is True
                if enforce_pairwise:
                    for existing_combo in valid_combinations:
                        pair_set.add(frozenset(set(combo) | set(existing_combo)))

        if len(valid_combinations) < GENES:
            raise ValueError(f"Unable to generate {GENES} unique codebook combinations with the constraints. "
                            f"Generated only {len(valid_combinations)} valid combinations.")

        best_codebook = None
        best_balance_score = float('inf')
        best_avg_expressions = {}
        
        num_iterations = min(10000, len(valid_combinations) ** 2)
        
        for iteration in range(num_iterations):
            shuffled_combinations = valid_combinations.copy()
            random.shuffle(shuffled_combinations)
            temp_codebook = shuffled_combinations[:GENES]
            
            number_expressions = {i: [] for i in range(1, n_numbers + 1)}
            
            for combo, expr in zip(temp_codebook, gene_expressions):
                for num in combo:
                    number_expressions[num].append(expr)
            
            avg_expressions = {}
            for num in range(1, n_numbers + 1):
                if number_expressions[num]:
                    avg_expressions[num] = np.mean(number_expressions[num])
                else:
                    avg_expressions[num] = 0
            
            balance_score = np.std(list(avg_expressions.values()))
            
            if balance_score < best_balance_score:
                best_balance_score = balance_score
                best_codebook = temp_codebook
                best_avg_expressions = avg_expressions

        print(f"Average normalized expression levels for each {readout_label} in the optimized codebook:")
        for i in range(1, n_numbers + 1):
            print(f"{readout_label} {i}: {best_avg_expressions[i]:.4f}")
        
        print(f"\nBalance score (standard deviation): {best_balance_score:.4f}")
        
        frequencies = np.zeros(n_numbers, dtype=int)
        for combo in best_codebook:
            for num in combo:
                frequencies[num - 1] += 1
        
        print(f"\nFrequencies of each {readout_label} in the codebook:")
        for i in range(1, n_numbers + 1):
            print(f"{readout_label} {i}: {frequencies[i - 1]}")

        codebook_df = pd.DataFrame(best_codebook, columns=[f'Read{i+1}' for i in range(READS)])
        return codebook_df

    def generate_binary_codes_with_exact_ones(HYBRIDIZATIONS, exact_ones):
        """Generate all binary codes with exactly exact_ones 1s.
        Mirrors the combinatorial logic where each barcode has exactly READS entries."""
        from itertools import combinations
        
        exact_ones = min(exact_ones, HYBRIDIZATIONS)
        all_codes = []
        for positions in combinations(range(HYBRIDIZATIONS), exact_ones):
            code = [0] * HYBRIDIZATIONS
            for pos in positions:
                code[pos] = 1
            all_codes.append(code)
        
        return all_codes
    
    def create_binary_codebook_per_channel(HYBRIDIZATIONS, CHANNELS, GENES, gene_expressions):
        """Generate per-channel binary codebooks with exactly READS ones, balanced within each channel.
        Mirrors combinatorial mode: pairwise uniqueness (union of active positions must be unique
        for every pair within a channel) and balancing are both applied.
        
        Returns:
            per_channel_barcodes, per_channel_assignments, genes_per_channel
        """
        # Generate all possible binary codes with exactly READS ones
        # (mirrors combinatorial mode where each barcode has exactly READS entries)
        # Use tuples so codes are hashable and membership checks are unambiguous.
        all_possible_codes = [tuple(c) for c in generate_binary_codes_with_exact_ones(HYBRIDIZATIONS, READS)]
        
        codes_per_channel_needed = (GENES + CHANNELS - 1) // CHANNELS  # round-up division
        
        if len(all_possible_codes) < codes_per_channel_needed:
            raise ValueError(f"Cannot generate {codes_per_channel_needed} unique binary barcodes per channel "
                           f"with {HYBRIDIZATIONS} imaging rounds and exactly {READS} ones per barcode. "
                           f"Maximum possible: {len(all_possible_codes)} (C({HYBRIDIZATIONS},{READS})). "
                           f"Either increase imaging rounds or reduce number of genes per channel.")
        
        use_expression_balancing = any(e != 1.0 for e in gene_expressions)
        
        # Distribute genes evenly across channels
        genes_per_channel = []
        gene_expressions_per_channel = []
        genes_assigned = 0
        
        for channel in range(CHANNELS):
            genes_in_this_channel = (GENES - genes_assigned + CHANNELS - channel - 1) // (CHANNELS - channel)
            genes_per_channel.append(genes_in_this_channel)
            
            channel_expressions = []
            for g_idx in range(genes_in_this_channel):
                global_gene_idx = genes_assigned + g_idx
                channel_expressions.append(gene_expressions[global_gene_idx] if global_gene_idx < len(gene_expressions) else 1.0)
            
            gene_expressions_per_channel.append(channel_expressions)
            genes_assigned += genes_in_this_channel
        
        # Helper: active bit positions of a code as a frozenset (for pairwise uniqueness)
        def active_positions(code):
            return frozenset(i for i, b in enumerate(code) if b == 1)
        
        # Helper: check and update pairwise uniqueness for a candidate code
        def pairwise_ok(code, selected_codes, pair_set):
            if not enforce_pairwise:
                return True
            code_pos = active_positions(code)
            for existing in selected_codes:
                if frozenset(code_pos | active_positions(existing)) in pair_set:
                    return False
            return True
        
        def register_pairs(code, selected_codes, pair_set):
            """Add all new pairs introduced by adding code to selected_codes (code already appended)."""
            if enforce_pairwise:
                code_pos = active_positions(code)
                for existing in selected_codes:  # includes code itself (self-pair)
                    pair_set.add(frozenset(code_pos | active_positions(existing)))
        
        # Generate balanced barcodes for each channel
        per_channel_barcodes = {}
        per_channel_assignments = {}
        
        def _run_selection_for_channel(num_genes_in_channel, channel_expressions, pairwise_enforced):
            """Run one selection attempt for a channel; returns selected_codes or None on failure.

            Uses its own pairwise check so it never mutates the outer enforce_pairwise variable.
            """
            shuffled_codes = all_possible_codes.copy()
            random.shuffle(shuffled_codes)

            selected_codes = []
            pair_set = set()

            def _pairwise_ok(code, sel_codes, p_set):
                if not pairwise_enforced:
                    return True
                code_pos = active_positions(code)
                seen_unions = set()
                for existing in sel_codes:
                    pair_union = frozenset(code_pos | active_positions(existing))
                    if pair_union in p_set:
                        return False
                    if pair_union in seen_unions:
                        return False
                    seen_unions.add(pair_union)
                return True

            def _register_pairs(code, sel_codes, p_set):
                if pairwise_enforced:
                    code_pos = active_positions(code)
                    for existing in sel_codes:
                        p_set.add(frozenset(code_pos | active_positions(existing)))

            if use_expression_balancing:
                bit_expressions = [[] for _ in range(HYBRIDIZATIONS)]

                for gene_pos in range(num_genes_in_channel):
                    gene_expr = channel_expressions[gene_pos] if gene_pos < len(channel_expressions) else 1.0
                    best_code = None
                    best_balance_score = float('inf')

                    for code in shuffled_codes:
                        if code in selected_codes:
                            continue
                        if not _pairwise_ok(code, selected_codes, pair_set):
                            continue

                        # Evaluate balance if this code were added
                        temp_bit_expressions = [list(be) for be in bit_expressions]
                        for bit_pos, bit in enumerate(code):
                            if bit == 1:
                                temp_bit_expressions[bit_pos].append(gene_expr)

                        avg_expressions = [np.mean(be) if be else 0.0 for be in temp_bit_expressions]
                        balance_score = np.std(avg_expressions)

                        if balance_score < best_balance_score:
                            best_balance_score = balance_score
                            best_code = code

                    if best_code is None:
                        return None  # failed to place this gene
                    selected_codes.append(best_code)
                    _register_pairs(best_code, selected_codes, pair_set)
                    for bit_pos, bit in enumerate(best_code):
                        if bit == 1:
                            bit_expressions[bit_pos].append(gene_expr)
            else:
                bit_frequencies = [0] * HYBRIDIZATIONS

                for gene_pos in range(num_genes_in_channel):
                    best_code = None
                    best_frequency_balance = float('inf')

                    for code in shuffled_codes:
                        if code in selected_codes:
                            continue
                        if not _pairwise_ok(code, selected_codes, pair_set):
                            continue

                        temp_frequencies = bit_frequencies.copy()
                        for bit_pos, bit in enumerate(code):
                            if bit == 1:
                                temp_frequencies[bit_pos] += 1

                        balance_score = np.std(temp_frequencies)

                        if balance_score < best_frequency_balance:
                            best_frequency_balance = balance_score
                            best_code = code

                    if best_code is None:
                        return None  # failed to place this gene
                    selected_codes.append(best_code)
                    _register_pairs(best_code, selected_codes, pair_set)
                    for bit_pos, bit in enumerate(best_code):
                        if bit == 1:
                            bit_frequencies[bit_pos] += 1

            return selected_codes

        MAX_RETRIES = 10

        for channel in range(CHANNELS):
            num_genes_in_channel = genes_per_channel[channel]
            channel_expressions = gene_expressions_per_channel[channel]

            if num_genes_in_channel == 0:
                continue

            selected_codes = None

            # Try up to MAX_RETRIES times with different random shuffles
            for attempt in range(MAX_RETRIES):
                selected_codes = _run_selection_for_channel(
                    num_genes_in_channel, channel_expressions, pairwise_enforced=enforce_pairwise
                )
                if selected_codes is not None and len(selected_codes) == num_genes_in_channel:
                    break
                selected_codes = None

            # If still failing and pairwise was enabled, retry with relaxed constraints
            if selected_codes is None and enforce_pairwise:
                print(f"  Warning: Could not satisfy pairwise uniqueness for channel {channel + 1} "
                      f"after {MAX_RETRIES} attempts. Retrying with relaxed pairwise constraints...")
                for attempt in range(MAX_RETRIES):
                    selected_codes = _run_selection_for_channel(
                        num_genes_in_channel, channel_expressions, pairwise_enforced=False
                    )
                    if selected_codes is not None and len(selected_codes) == num_genes_in_channel:
                        print(f"  Warning: Channel {channel + 1} codebook generated without pairwise "
                              f"uniqueness enforcement.")
                        break
                    selected_codes = None

            if selected_codes is None or len(selected_codes) != num_genes_in_channel:
                assigned = len(selected_codes) if selected_codes else 0
                raise ValueError(
                    f"Binary barcode generation failed for channel {channel + 1}: "
                    f"could only assign {assigned} barcodes to {num_genes_in_channel} genes. "
                    f"Not enough distinct binary codes exist with {HYBRIDIZATIONS} imaging rounds "
                    f"and exactly {READS} ones to satisfy all constraints. "
                    f"Try increasing the number of imaging rounds or reducing the number of genes."
                )

            # Convert back to lists for downstream compatibility
            per_channel_barcodes[channel] = [list(c) for c in selected_codes]
            per_channel_assignments[channel] = {'genes_count': num_genes_in_channel}

            print(f"\nChannel {channel + 1} - {num_genes_in_channel} genes:")
            for pos in range(HYBRIDIZATIONS):
                freq = sum(1 for code in selected_codes if code[pos] == 1)
                print(f"  Imaging round {pos + 1}: {freq} genes")
        
        return per_channel_barcodes, per_channel_assignments, genes_per_channel

    # Choose which method to use based on barcode_type parameter
    if barcode_type == "binary":
        print("Using binary barcode generation (per-channel with max 4 ones)")
        per_channel_barcodes, per_channel_assignments, genes_per_channel = \
            create_binary_codebook_per_channel(HYBRIDIZATIONS, CHANNELS, GENES, gene_expressions)
        
        # Flatten the assignment for easier use later
        all_channel_assignments = []
        all_binary_codes = []
        
        for channel in range(CHANNELS):
            if channel in per_channel_barcodes:
                for code in per_channel_barcodes[channel]:
                    all_channel_assignments.append(channel + 1)  # 1-indexed
                    all_binary_codes.append(code)
        
        binary_codes = all_binary_codes
        channel_assignments = all_channel_assignments
        is_binary_mode = True
    else:
        # Original combinatorial barcode method
        if optimize_balance:
            # Original behavior: use expression-based optimization if available
            use_expression_optimization = any(e != 1.0 for e in gene_expressions)
            if use_expression_optimization:
                codebook_df = create_balanced_codebook(READS, n_numbers, GENES, gene_expressions)
            else:
                codebook_df = create_max_codebook(READS, n_numbers, GENES)
        else:
            # Simple method: maximize valid barcodes without iterative optimization
            print("Using simple codebook generation (no frequency optimization)")
            codebook_df = create_simple_codebook(READS, n_numbers, GENES)
        is_binary_mode = False

    # ===================== BINARY MODE =====================
    if is_binary_mode:
        print("\n=== Binary Barcode Processing (Per-Channel) ===")
        
        # Reconstruct gene assignments to channels based on per_channel_barcodes
        gene_to_channel_barcode = {}  # gene_idx -> (channel, barcode_index_in_channel)
        gene_idx = 0
        
        for channel in range(CHANNELS):
            if channel in per_channel_barcodes:
                for barcode_idx, barcode in enumerate(per_channel_barcodes[channel]):
                    if gene_idx < GENES:
                        gene_to_channel_barcode[gene_idx] = (channel + 1, barcode_idx, barcode)
                        gene_idx += 1
        
        # Generate pattern mapping for binary mode
        # For each (hybridization, channel) pair, assign a sequence
        pattern_data = []
        sequence_idx = 0
        
        for hyb_round in range(1, HYBRIDIZATIONS + 1):
            for channel in range(1, CHANNELS + 1):
                if sequence_idx < len(strings):
                    pattern_data.append({
                        'Hybridization': hyb_round,
                        'Channel': channel,
                        'SequenceIndex': sequence_idx,
                        'Sequence': strings[sequence_idx]
                    })
                    sequence_idx += 1
        
        # Create binary codebook output dataframe
        binary_output_columns = ['Gene', 'Channel'] + [f'H{i+1}' for i in range(HYBRIDIZATIONS)]
        binary_codebook_output = []
        
        binary_sequences_for_probes = []  # For storing sequences used in barcodes
        
        for gene_idx, gene in enumerate(genes):
            if gene_idx in gene_to_channel_barcode:
                channel, barcode_idx, binary_code = gene_to_channel_barcode[gene_idx]
                
                row = [gene, channel]
                row.extend(binary_code)
                binary_codebook_output.append(row)
                
                # Collect sequences for this gene based on binary code and channel
                gene_sequences = []
                for hyb_idx, bit in enumerate(binary_code):
                    if bit == 1:
                        # Find sequence for this (hyb_idx+1, channel) pair
                        seq_idx = None
                        for pattern in pattern_data:
                            if pattern['Hybridization'] == hyb_idx + 1 and pattern['Channel'] == channel:
                                seq_idx = pattern['SequenceIndex']
                                break
                        
                        if seq_idx is not None and seq_idx < len(strings):
                            gene_sequences.append(strings[seq_idx])
                
                binary_sequences_for_probes.append(gene_sequences)
        
        binary_codebook_df = pd.DataFrame(binary_codebook_output, columns=binary_output_columns)
        
        # For binary mode, create a simplified pattern DataFrame
        pattern_df = pd.DataFrame(pattern_data)
        if MODE == "RAEFISH":
            pattern_df = pattern_df[['Hybridization', 'Channel', 'Sequence']].copy()
            pattern_df.columns = ['Hybridization', 'Channel', 'Adapter_sequence']
        else:
            pattern_df = pattern_df[['Hybridization', 'Channel', 'Sequence']].copy()
            pattern_df.columns = ['Hybridization', 'Channel', 'Readout_sequence']
        
        # Handle output for binary MERFISH
        if MODE == "MERFISH":
            print(f"Generating MERFISH binary codebook with {len(genes)} genes, {CHANNELS} channels, and {HYBRIDIZATIONS} imaging rounds")
            
            # For binary MERFISH, encode the binary position information in the probes
            five_sec_seq = []
            three_sec_seq = []
            
            for gene_idx, gene_sequences in enumerate(binary_sequences_for_probes):
                # Each gene must have exactly READS sequences (one per '1' bit in barcode)
                if len(gene_sequences) != READS:
                    raise ValueError(f"Gene '{genes[gene_idx]}' has {len(gene_sequences)} active imaging rounds but READS={READS}. "
                                     f"Binary barcodes must have exactly {READS} ones.")
                
                s_values = gene_sequences.copy()
                random.shuffle(s_values)
                
                def random_letters(n):
                    nucleotides = ['A', 'T', 'C', 'G']
                    result = []
                    last_nt = None
                    for _ in range(n):
                        choices = [nt for nt in nucleotides if nt != last_nt]
                        new_nt = random.choice(choices)
                        result.append(new_nt)
                        last_nt = new_nt
                    return ''.join(result)
                
                five_sec_seq_str = f"{FIVE_PRIMER}{s_values[0]}{random_letters(1)}{s_values[1]}{random_letters(2)}"
                three_sec_seq_str = f"{random_letters(2)}{s_values[2]}{random_letters(1)}{s_values[3]}{THREE_PRIMER}"
                
                five_sec_seq.append(five_sec_seq_str)
                three_sec_seq.append(three_sec_seq_str)
            
            if reverse_complement:
                five_sec_seq = [str(Seq(seq).reverse_complement()) for seq in five_sec_seq]
                three_sec_seq = [str(Seq(seq).reverse_complement()) for seq in three_sec_seq]
                five_sec_seq, three_sec_seq = three_sec_seq, five_sec_seq
            
            sec_sequences_df = pd.DataFrame({
                '5_sec_sequences': five_sec_seq,
                '3_sec_sequences': three_sec_seq
            })
            
            to_order_idt_df = _build_merfish_to_order_idt(mouse_seq_file, pattern_df)
            with pd.ExcelWriter(OUTPUT_FILE) as writer:
                pattern_df.to_excel(writer, sheet_name='Sequence_Pattern', index=False)
                binary_codebook_df.to_excel(writer, sheet_name='Gene_Codebook', index=False)
                sec_sequences_df.to_excel(writer, sheet_name='Encoding_probe_arms', index=False)
                to_order_idt_df.to_excel(writer, sheet_name='To_order_IDT', index=False)
        
        elif MODE == "RAEFISH":
            print(f"Generating RAEFISH binary codebook with {len(genes)} genes, {CHANNELS} channels, and {HYBRIDIZATIONS} imaging rounds")
            
            # For binary RAEFISH, create adapter probes
            adapter_probes_data = []
            for hyb_round in range(1, HYBRIDIZATIONS + 1):
                for channel in range(1, CHANNELS + 1):
                    # Find sequence for this pair
                    seq_to_use = None
                    for idx, pattern in enumerate(pattern_data):
                        if pattern['Hybridization'] == hyb_round and pattern['Channel'] == channel:
                            seq_to_use = pattern['Sequence']
                            break
                    
                    if seq_to_use:
                        original_adapter = seq_to_use  # Already the adapter sequence
                        landing_sequence = raefish_readout_probes.get(channel, '')
                        adapter_probe = f"{original_adapter}{landing_sequence}{landing_sequence}"
                        adapter_probes_data.append({
                            'Hybridization': hyb_round,
                            'Channel': channel,
                            'Adapter_sequence': seq_to_use,
                            'Adapter_probe': adapter_probe
                        })
            
            pattern_df_binary = pd.DataFrame(adapter_probes_data)
            
            # Create encoding probes for binary RAEFISH
            five_sec_seq = []
            three_sec_seq = []
            
            def random_letters(n):
                nucleotides = ['A', 'T', 'C', 'G']
                result = []
                last_nt = None
                for _ in range(n):
                    choices = [nt for nt in nucleotides if nt != last_nt]
                    new_nt = random.choice(choices)
                    result.append(new_nt)
                    last_nt = new_nt
                return ''.join(result)
            
            for gene_idx, gene_sequences in enumerate(binary_sequences_for_probes):
                # Each gene must have exactly READS sequences (one per '1' bit in barcode)
                if len(gene_sequences) != READS:
                    raise ValueError(f"Gene '{genes[gene_idx]}' has {len(gene_sequences)} active imaging rounds but READS={READS}. "
                                     f"Binary barcodes must have exactly {READS} ones.")
                
                s_values = gene_sequences.copy()
                random.shuffle(s_values)
                
                five_sec_seq_str = f"{FIVE_PRIMER}{s_values[0]}{random_letters(1)}{s_values[1]}{random_letters(2)}"
                three_sec_seq_str = f"{random_letters(2)}{s_values[2]}{random_letters(1)}{s_values[3]}{THREE_PRIMER}"
                
                five_sec_seq.append(five_sec_seq_str)
                three_sec_seq.append(three_sec_seq_str)
            
            # Create encoding probes for RAEFISH
            encoding_probe1_list = []
            encoding_probe2_list = []
            padlock_1_list = []
            padlock_2_list = []
            priming_1_list = []
            priming_2_list = []
            
            for i, gene in enumerate(genes):
                if gene in raefish_gene_mapping:
                    target_1 = raefish_gene_mapping[gene]['Target_1']
                    target_2 = raefish_gene_mapping[gene]['Target_2']
                    padlock_1 = raefish_gene_mapping[gene]['Padlock_1']
                    padlock_2 = raefish_gene_mapping[gene]['Padlock_2']
                    priming_1 = raefish_gene_mapping[gene]['Priming_1']
                    priming_2 = raefish_gene_mapping[gene]['Priming_2']
                    
                    # Convert to str and treat NaN (empty Excel cells) as missing
                    import math
                    def _is_missing(v):
                        return v is None or (isinstance(v, float) and math.isnan(v))
                    
                    encoding_probe1 = f"{five_sec_seq[i]}{target_1}{three_sec_seq[i]}" if not _is_missing(target_1) else None
                    encoding_probe2 = f"{five_sec_seq[i]}{target_2}{three_sec_seq[i]}" if not _is_missing(target_2) else None
                    
                    encoding_probe1_list.append(encoding_probe1)
                    encoding_probe2_list.append(encoding_probe2)
                    padlock_1_list.append(None if _is_missing(padlock_1) else padlock_1)
                    padlock_2_list.append(None if _is_missing(padlock_2) else padlock_2)
                    priming_1_list.append(None if _is_missing(priming_1) else priming_1)
                    priming_2_list.append(None if _is_missing(priming_2) else priming_2)
                else:
                    print(f"Warning: Gene '{gene}' not found in RAEFISH_sequences.")
                    encoding_probe1_list.append(None)
                    encoding_probe2_list.append(None)
                    padlock_1_list.append(None)
                    padlock_2_list.append(None)
                    priming_1_list.append(None)
                    priming_2_list.append(None)
            
            encoding_probes_rc_list1 = [str(Seq(probe).reverse_complement()) if probe else None for probe in encoding_probe1_list]
            encoding_probes_rc_list2 = [str(Seq(probe).reverse_complement()) if probe else None for probe in encoding_probe2_list]

            primary_probes_df = pd.DataFrame({
                'Gene': genes,
                'Encoding_probe1': encoding_probe1_list,
                'Encoding_probe2': encoding_probe2_list,
                'Padlock_1': padlock_1_list,
                'Padlock_2': padlock_2_list,
                'Priming_1': priming_1_list,
                'Priming_2': priming_2_list
            })
            to_order_twist_df = _build_raefish_to_order_twist(
                genes, encoding_probes_rc_list1, encoding_probes_rc_list2,
                padlock_1_list, padlock_2_list, priming_1_list, priming_2_list)
            to_order_idt_df = _build_raefish_to_order_idt(mouse_seq_file, pattern_df_binary)

            with pd.ExcelWriter(OUTPUT_FILE) as writer:
                pattern_df_binary.to_excel(writer, sheet_name='Sequence_Pattern', index=False)
                binary_codebook_df.to_excel(writer, sheet_name='Gene_Codebook', index=False)
                primary_probes_df.to_excel(writer, sheet_name='Primary_Probes', index=False)
                to_order_twist_df.to_excel(writer, sheet_name='To_order_Twist', index=False)
                to_order_idt_df.to_excel(writer, sheet_name='To_order_IDT', index=False)
        
        print(f"\nBinary codebook generation complete for {MODE} mode")
        print(f"Output saved to {OUTPUT_FILE}")
        return binary_codebook_df, pattern_df
    
    # ===================== COMBINATORIAL MODE =====================
    # Step 3: Generate Pattern DataFrame
    col1, col2, col3 = [], [], []
    readout_value = 1

    for j in range(1, HYBRIDIZATIONS + 1):
        for k in range(1, CHANNELS + 1):
            col1.append(j)
            col2.append(k)
            col3.append(readout_value)
            readout_value += 1

    if MODE == "RAEFISH":
        pattern_df = pd.DataFrame({
            'Adapter_number': range(1, min_sec_seq + 1),
            'Adapter_sequence': strings[:min_sec_seq],
            'Hybridization': col1,
            'Channel': col2,
            'Readout': col3
        })
    else:  # MERFISH
        pattern_df = pd.DataFrame({
            'Readout_number': range(1, min_sec_seq + 1),
            'Readout_sequence': strings[:min_sec_seq],
            'Hybridization': col1,
            'Channel': col2,
            'Readout': col3
        })

    # Step 4: Place Strings in Table Based on Codebook and Pattern Columns
    placed_strings = pd.DataFrame(index=range(len(codebook_df)), columns=[f'Pos{col}' for col in range(1, READS + 1)])

    # Determine the sequence column name based on mode
    if MODE == "RAEFISH":
        sequence_column = 'Adapter_sequence'
    else:  # MERFISH
        sequence_column = 'Readout_sequence'

    for row_codebook in range(len(codebook_df)):
        for read_num in range(1, READS + 1):
            readout_num = codebook_df.iloc[row_codebook, read_num - 1]
            matching_pattern_row = pattern_df[pattern_df['Readout'] == readout_num]
            if not matching_pattern_row.empty:
                string_value = matching_pattern_row.iloc[0][sequence_column]
                placed_strings.iloc[row_codebook, read_num - 1] = string_value

    # Step 5: Create Secondary Sequences DataFrame with Shuffling
    def random_letters(n):
        """Generate a random sequence of 'n' nucleotides with no consecutive repeats."""
        nucleotides = ['A', 'T', 'C', 'G']
        result = []
        last_nt = None
        for _ in range(n):
            choices = [nt for nt in nucleotides if nt != last_nt]
            new_nt = random.choice(choices)
            result.append(new_nt)
            last_nt = new_nt
        return ''.join(result)
 
    five_sec_seq, three_sec_seq = [], []

    for i in range(len(codebook_df)):
        s_values = placed_strings.iloc[i].tolist()
        random.shuffle(s_values)

        five_sec_seq_str = f"{FIVE_PRIMER}{s_values[0]}{random_letters(1)}{s_values[1]}{random_letters(2)}"
        three_sec_seq_str = f"{random_letters(2)}{s_values[2]}{random_letters(1)}{s_values[3]}{THREE_PRIMER}"

        five_sec_seq.append(five_sec_seq_str)
        three_sec_seq.append(three_sec_seq_str)

    # Step 6: Handle mode-specific processing
    if MODE == "MERFISH":
        # Conditionally apply reverse complement based on reverse_complement parameter
        if reverse_complement:
            five_sec_seq = [str(Seq(seq).reverse_complement()) for seq in five_sec_seq]
            three_sec_seq = [str(Seq(seq).reverse_complement()) for seq in three_sec_seq]
            five_sec_seq, three_sec_seq = three_sec_seq, five_sec_seq

        sec_sequences_df = pd.DataFrame({
            '5_sec_sequences': five_sec_seq,
            '3_sec_sequences': three_sec_seq
        })

        # Create output with gene names
        codebook_output = codebook_df.copy()
        codebook_output.insert(0, 'Gene', genes)

        # Write to Excel
        to_order_idt_df = _build_merfish_to_order_idt(mouse_seq_file, pattern_df)
        with pd.ExcelWriter(OUTPUT_FILE) as writer:
            pattern_df.to_excel(writer, sheet_name='Sequence_Pattern', index=False)
            codebook_output.to_excel(writer, sheet_name='Gene_Codebook', index=False)
            sec_sequences_df.to_excel(writer, sheet_name='Encoding_probe_arms', index=False)
            to_order_idt_df.to_excel(writer, sheet_name='To_order_IDT', index=False)

    elif MODE == "RAEFISH":
        # Create adapter probes for RAEFISH
        adapter_probes_list = []
        
        for i in range(len(pattern_df)):
            # Get the original adapter sequence (before reverse complement)
            original_adapter = raefish_adapter_sequences[i]
            channel = pattern_df.iloc[i]['Channel']
            
            # Get landing sequence for this channel
            landing_sequence = raefish_readout_probes[channel]
            
            # Create adapter probe: original_adapter + landing_sequence landing_sequence
            adapter_probe = f"{original_adapter}{landing_sequence}{landing_sequence}"
            adapter_probes_list.append(adapter_probe)
        
        # Update pattern_df with adapter probes
        pattern_df_raefish = pattern_df.copy()
        pattern_df_raefish.insert(2, 'Adapter_probe', adapter_probes_list)

        # Create encoding probes and additional sequences for RAEFISH
        encoding_probe1_list = []
        encoding_probe2_list = []
        padlock_1_list = []
        padlock_2_list = []
        priming_1_list = []
        priming_2_list = []

        import math
        def _is_missing(v):
            return v is None or (isinstance(v, float) and math.isnan(v))
        
        for i, gene in enumerate(genes):
            if gene in raefish_gene_mapping:
                target_1 = raefish_gene_mapping[gene]['Target_1']
                target_2 = raefish_gene_mapping[gene]['Target_2']
                padlock_1 = raefish_gene_mapping[gene]['Padlock_1']
                padlock_2 = raefish_gene_mapping[gene]['Padlock_2']
                priming_1 = raefish_gene_mapping[gene]['Priming_1']
                priming_2 = raefish_gene_mapping[gene]['Priming_2']
                
                # Convert to str and treat NaN (empty Excel cells) as missing
                encoding_probe1 = f"{five_sec_seq[i]}{target_1}{three_sec_seq[i]}" if not _is_missing(target_1) else None
                encoding_probe2 = f"{five_sec_seq[i]}{target_2}{three_sec_seq[i]}" if not _is_missing(target_2) else None
                
                encoding_probe1_list.append(encoding_probe1)
                encoding_probe2_list.append(encoding_probe2)
                padlock_1_list.append(None if _is_missing(padlock_1) else padlock_1)
                padlock_2_list.append(None if _is_missing(padlock_2) else padlock_2)
                priming_1_list.append(None if _is_missing(priming_1) else priming_1)
                priming_2_list.append(None if _is_missing(priming_2) else priming_2)
            else:
                print(f"Warning: Gene '{gene}' not found in RAEFISH_sequences.")
                encoding_probe1_list.append(None)
                encoding_probe2_list.append(None)
                padlock_1_list.append(None)
                padlock_2_list.append(None)
                priming_1_list.append(None)
                priming_2_list.append(None)

        # Create reverse complement dataframe for ordering
        encoding_probes_rc_list1 = [str(Seq(probe).reverse_complement()) if probe else None for probe in encoding_probe1_list]
        encoding_probes_rc_list2 = [str(Seq(probe).reverse_complement()) if probe else None for probe in encoding_probe2_list]

        primary_probes_df = pd.DataFrame({
            'Gene': genes,
            'Encoding_probe1': encoding_probe1_list,
            'Encoding_probe2': encoding_probe2_list,
            'Padlock_1': padlock_1_list,
            'Padlock_2': padlock_2_list,
            'Priming_1': priming_1_list,
            'Priming_2': priming_2_list
        })
        to_order_twist_df = _build_raefish_to_order_twist(
            genes, encoding_probes_rc_list1, encoding_probes_rc_list2,
            padlock_1_list, padlock_2_list, priming_1_list, priming_2_list)
        to_order_idt_df = _build_raefish_to_order_idt(mouse_seq_file, pattern_df_raefish)

        # Create output with gene names
        codebook_output = codebook_df.copy()
        codebook_output.insert(0, 'Gene', genes)

        # Write to Excel
        with pd.ExcelWriter(OUTPUT_FILE) as writer:
            pattern_df_raefish.to_excel(writer, sheet_name='Sequence_Pattern', index=False)
            codebook_output.to_excel(writer, sheet_name='Gene_Codebook', index=False)
            primary_probes_df.to_excel(writer, sheet_name='Primary_Probes', index=False)
            to_order_twist_df.to_excel(writer, sheet_name='To_order_Twist', index=False)
            to_order_idt_df.to_excel(writer, sheet_name='To_order_IDT', index=False)

    print(f"\nCodebook generation complete for {MODE} mode")
    print(f"Output saved to {OUTPUT_FILE}")
    return codebook_df, pattern_df


# IMPORTANT: Only run the codebook generation if this script is run directly, not when imported
if __name__ == "__main__":
    # Validate mode
    if MODE not in ["MERFISH", "RAEFISH"]:
        raise ValueError(f"MODE must be 'MERFISH' or 'RAEFISH', not '{MODE}'")

    # Load primers and readout sequences
    MOUSE_SEQUENCES_FILE = get_mouse_sequences_file()

    # Load gene list if using MERFISH or RAEFISH
    gene_expressions = None
    genes = None
    raefish_adapter_sequences = None
    raefish_readout_probes = None
    raefish_gene_mapping = None

    if MODE == "MERFISH":
        if GENE_LIST_FILE is None:
            raise ValueError("GENE_LIST_FILE must be specified for MERFISH mode")
        FIVE_PRIMER, THREE_PRIMER = load_primers_from_file(MOUSE_SEQUENCES_FILE, primer_type='MERFISH')
        genes, gene_expressions = load_gene_list(GENE_LIST_FILE)
        GENES = len(genes)
        print(f"Loaded {GENES} genes from {GENE_LIST_FILE}")
        print(f"Using {'expression-based' if any(e != 1.0 for e in gene_expressions) else 'standard'} optimization for codebook generation")
        READOUT_SEQUENCES_LIST = load_readout_sequences_from_file(MOUSE_SEQUENCES_FILE, sheet_name='MERFISH_readout_sequences')
        codebook_df, pattern_df = generate_hybridization_codebook(MODE, GENES, READS, HYBRIDIZATIONS, CHANNELS, READOUT_SEQUENCES_LIST, genes, gene_expressions, FIVE_PRIMER, THREE_PRIMER, OUTPUT_FILE, reverse_complement=False)
        
    elif MODE == "RAEFISH":
        if GENE_LIST_FILE is None:
            raise ValueError("GENE_LIST_FILE must be specified for RAEFISH mode")
        FIVE_PRIMER, THREE_PRIMER = load_primers_from_file(MOUSE_SEQUENCES_FILE, primer_type='RAEFISH')
        genes, gene_expressions = load_gene_list(GENE_LIST_FILE)
        GENES = len(genes)
        print(f"Loaded {GENES} genes from {GENE_LIST_FILE}")
        print(f"Using {'expression-based' if any(e != 1.0 for e in gene_expressions) else 'standard'} optimization for codebook generation")
        raefish_adapter_sequences = load_raefish_adapter_sequences(MOUSE_SEQUENCES_FILE)
        raefish_readout_probes = load_raefish_readout_probes(MOUSE_SEQUENCES_FILE)
        raefish_gene_mapping = load_raefish_sequences(MOUSE_SEQUENCES_FILE)
        # Reverse complement adapter sequences for RAEFISH
        READOUT_SEQUENCES_LIST = [str(Seq(seq).reverse_complement()) for seq in raefish_adapter_sequences]
        codebook_df, pattern_df = generate_hybridization_codebook(MODE, GENES, READS, HYBRIDIZATIONS, CHANNELS, READOUT_SEQUENCES_LIST, genes, gene_expressions, FIVE_PRIMER, THREE_PRIMER, OUTPUT_FILE, raefish_adapter_sequences, raefish_readout_probes, raefish_gene_mapping, readout_label="Adapter")
