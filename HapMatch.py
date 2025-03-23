#!/usr/bin/python3
import streamlit as st
import pandas as pd
from io import StringIO

# Sets the session stage in the overall process
if "stage" not in st.session_state:
    st.session_state.stage = 0
def set_stage(stage):
    st.session_state.stage = stage

st.markdown("## MtDNA Haplogroup Matching App :dna:")
st.markdown("### Jacob Daigle")
st.markdown("##### Lund University | March 2025")
st.markdown("***")
st.sidebar.markdown("""
                 # Notation :memo:\n
                 "+" both individuals have the mutation.\n
                 "-" neither individual has the mutation.\n
                 "/" only the second (reference) individual has the mutation.\n
                 "\\\\" only the first (input) individual has the mutation.\n
                "d" at least one individual has a deletion.\n
                "?" at least one individual has an ambiguous nucleotide base.
                 """)
st.sidebar.markdown("# Options :gear:")
st.session_state["show_map"] = st.sidebar.checkbox(label = "Show map")
READ_LIMIT = st.sidebar.number_input(
    "Reference Read Limit",
    min_value=0,
    value=0,
    step=1,
    help="Set a max number of reference sequences to parse to save on computation time. Set to 0 for no limit."
    )
st.session_state.comparison_mode = st.radio(
    "-> Choose what type of comparison you'd like to make.",
    ("User vs. Ancient", "User vs. User", "Ancient vs. Ancient"),
    horizontal=True,
    on_change=set_stage,
    args=(0,)
    )
st.write("")
st.write("-> Upload a user file, enter the user's haplogroup, and hit run!")
input_box = st.empty()
error_text = st.empty()
table = st.empty()

# Nucleotide compare: compares nucleotide bases including ambiguities
symbols = {
    "A": ("A"),
    "T": ("T"),
    "C": ("C"),
    "G": ("G"),
    "U": ("U", "T"),
    "M": ("A", "C"),
    "R": ("A", "G"),
    "W": ("A", "T"),
    "S": ("C", "G"),
    "Y": ("C", "T"),
    "K": ("G", "T"),
    "H": ("A", "C", "T"),
    "D": ("A", "G", "T"),
    "V": ("A", "C", "G"),
    "B": ("C", "G", "T"),
    "N": ("A", "T", "C", "G"),
    "X": ("A", "T", "C", "G"),
    "-": ("-"), "--": ("--"),
    "I": ("A")
}
ambiguous = ("M", "R", "W", "S", "Y", "K", "H", "D", "V", "B", "N", "X")
deletion = ("-", "--")
def nc(base1, base2):
    base1 = base1.upper()
    base2 = base2.upper()
    if base1 not in symbols or base2 not in symbols:  # Catches any more non-standard symbols
        return False
    return (base1 in symbols[base2]) or (base2 in symbols[base1])


# Compares user file to specific reference individual
@st.cache_data
def compareTo(usr, seq, snps: dict, mode: str) -> tuple:
    not_found = set()
    found = list()
    # haplotypes.items looks like [(123, (g, c123g)), ...] = [(pos, (base, id)), ...]
    for pos, mutation in snps.items():
        base = mutation[0]
        id = mutation[1]
        if mode == "ua":
            if pos not in usr:
                not_found.add(f"<{id}>")
                continue
            else:
                usr_base = usr[pos]
                ref_base = seq[pos-1]
        elif mode == "aa":
            usr_base = usr[pos-1]
            ref_base = seq[pos-1]
        elif mode == "uu":
            if pos not in usr or pos not in seq:
                not_found.add(f"<{id}>")
                continue
            else:
                usr_base = usr[pos]
                ref_base = seq[pos]

        ref_usr_matching = (nc(ref_base, base), nc(usr_base, base))
        statement = ""
        match ref_usr_matching:
            case (False, False):
                statement = f"-<{id}>"  # If neither reference nor user have mutation, it is -
            case (True, True):
                statement = f"+<{id}>"  # If both have the mutation
            case (True, False):
                statement = f"\\<{id}>"  # If only reference has the mutation
            case (False, True):
                statement = f"/<{id}>"  # If only user has the mutation
        if ref_base in ambiguous or usr_base in ambiguous:
            statement += "?"
        if ref_base in deletion or usr_base in deletion:
            statement += "d"
        found.append(statement)
    return (";".join(found), not_found)


# Compares user file to all reference individuals and scores based on differences
@st.cache_data
def getClosest(refs: dict, usr, haplos: dict, mode: str) -> dict:
    mode = {"User vs. Ancient": "ua", "User vs. User": "uu", "Ancient vs. Ancient": "aa"}[mode]
    if mode == "aa":
        usr = refs[usr]
    # Get number of differences (score) for each reference
    similarity = dict()
    for id, seq in refs.items():
        differences = 0
        for pos in haplos:
            if mode == "ua":
                if pos in usr and not nc(seq[pos - 1], usr[pos]):
                    differences += 1
            elif mode == "aa":
                differences += int(not nc(seq[pos-1], usr[pos-1]))
            elif mode == "uu":
                if pos in usr and pos in seq and not nc(usr[pos], seq[pos]):  # In uu case, id is User 2 and seq is {pos: base} dict
                    print(usr[pos], seq[pos])
                    differences += 1
        similarity[id] = [differences]  # List because mutations will be appended
    
    # Go through each id and compare to get present/missing mutations
    for id in similarity:
        missing = set()
        results = compareTo(usr, refs[id], haplos, mode)  # ("+<mutation1>;-<mutation2>;...", "<mutationX>;<mutationY>;...")
        if results[1]:
            missing |= results[1]
        if missing:
            pass
        similarity[id].append(results[0])

    # Convert dictionary to pandas dataframe
    df = pd.DataFrame.from_dict(
        data=similarity,
        orient="index",
        columns=["Genetic Distance", "Mutations"]
        )
    df.index.name = "Name"
    df.sort_values("Genetic Distance", inplace=True)

    return df


# Parses mutations file into haplo[pos] = (base, name)
@st.cache_data
def get_haplo(path: str, ids: list) -> dict:
    try:
        haplo = dict()
        with open(path) as f:
            for line in f:
                line = line.strip().split("\t")
                if line[3] in ids:
                    haplo[int(line[1])] = (line[2].casefold(), line[0])
        return haplo
    
    except FileNotFoundError:
        print(f"Error: file {path} not found.")


# Parses User file to get genome[position] = genotype
@st.cache_data
def parse_user(path="", file=None) -> dict:
    def parse(f):
        genome = dict()
        header_parsed = False
        for line in f:
                line = line.casefold().strip().split("\t")
                if line[0].startswith("#"):  # Parse the header
                    header_parsed = True
                    if "chromosome" in line:
                        chromosome_index = line.index("chromosome")
                    elif "chrom" in line:
                        chromosome_index = line.index("chrom")
                    else:
                        error_text.write("Parsing error: no chromosome/chrom column found in user file")
                        break

                    if "position" in line:
                        position_index = line.index("position")
                    elif "pos" in line:
                        position_index = line.index("pos")
                    else:
                        error_text.write("Parsing error: no position/pos column found in user file.")
                        break

                    if "genotype" in line:
                        genotype_index = line.index("genotype")
                    elif "allele1" in line:
                        genotype_index = line.index("allele1")
                    else:
                        error_text.write("Parsing error: no genotype/allele1 column found in user file.")
                        break
                elif header_parsed:
                    if line[chromosome_index] in ["mt", "0", "26"]:  # Only mitochondrial mutations
                        genome[int(line[position_index])] = line[genotype_index].strip()
                else:
                    error_text.write("Parsing error: no header row found in user file.")
                    break
        else:
            return genome
        return None

    if path:
        try:
            with open(path) as f:
                return parse(f)
        except FileNotFoundError:
            print(f"{path} not found.")
    elif file:
        return parse(StringIO(file.getvalue().decode("utf-8")))


# Parses reference file into parsed[id] = sequence
@st.cache_data
def parse_ref(path: str, limit=None) -> dict:
    try:
        with open(path) as f:
            total_seqs = f.read().count(">")
        if not limit:
            limit = total_seqs
        with open(path) as f:
            parsed = dict()
            line = f.readline()
            seq = ""
            header = ""
            count = 0
            while line:
                if line.startswith(">"):
                    count += 1
                    # updateProgress(count, limit)  # Progress bar
                    if count >= limit:  # Limiting number of reads to save computation time
                        break
                    elif header:  # header is empty on first entry only
                        parsed[header] = seq
                        seq = ""
                    header = line.strip()[1:]  # Removes >
                else:  # If not header, append section to seq until next header is reached
                    seq += line.strip().casefold()
                line = f.readline()

                if not line:  # This catches the last sequence in the file
                    parsed[header] = seq

        return parsed

    except FileNotFoundError:
        print(f"Error: file {path} not found.")


def show_table(display_box):
        return display_box.dataframe(st.session_state.results, on_select="rerun", selection_mode="single-row")


def get_metadata(names):
    global metadata_header
    if type(names) != list:
        names = [names]
    results = []
    with open("data/mtdb_metadata.tsv") as f:
            for line in f:
                if line.startswith("#"):
                    metadata_header = line.strip()[1:].split("\t")
                for name in names:
                    if name in line:
                        results.append(line.strip().split("\t"))
    if len(results) == 1:
        return results[0]
    else:
        return results


def main():
    match st.session_state.comparison_mode:
        case "User vs. Ancient":
            upload_box, haplo_box = input_box.columns(2)
            uploads = (upload_box,)
        case "User vs. User":
            upload_box1, upload_box2, haplo_box = input_box.columns(3)
            uploads = (upload_box1, upload_box2)
        case "Ancient vs. Ancient":
            ancient_user_select, haplo_display, run_button = input_box.columns([.7, .2, .1], vertical_alignment="bottom")
            haplo_box = None
            uploads = None

    reference_path = "data/mtdb.fasta"  # Fasta of ancient individuals
    references = parse_ref(reference_path, READ_LIMIT)

    users = []
    if uploads:
        user_files = []
        for i, box in enumerate(uploads):
            user_files.append(box.file_uploader(
                "Upload a user file",
                help="File should have headers/metadata marked by # and subsequent lines should be mutation statuses e.g. id | position | genotype",
                key=f"upload_box_{i}",
                on_change=set_stage,
                args = [0]
            ))
            # user_path = "Test1.tsv"
        for ufile in user_files:
            if ufile is not None:
                users.append(parse_user(file=ufile))
            # else:
            #     user = parse_user(path=user_path)
    else:
        users.append(
            ancient_user_select.selectbox(
                "Select an ancient individual",
                options=references.keys()
                ))
    
    if haplo_box:
        with haplo_box.container():
            main_group = st.text_input(  # Gets the user-input haplogroup
                "Haplogroup",
                help="A collection of haplotypes with common ancestry. E.g. X2c2, K1a11, N9a10b"
                )
            run_button = st.empty()
    else:
        main_group = get_metadata(users[0])[6]
        haplo_display.text_input(
            "Haplogroup",
            placeholder=main_group,
            disabled=True
            )

    subgroups = [main_group[:len(main_group)-i] for i in range(len(main_group))]  # E.g. K1a4b -> K1a4b, K1a4, K1a, K1, K
    haplo_path = "data/Mutation.txt"  # Haplogroup-SNP association file
    haplogroup = get_haplo(haplo_path, subgroups)

    run_button.button(
        label="Run",
        on_click=set_stage,
        args = [1],
        use_container_width=True,
        disabled=(None in users) or (not users)
        )

    if st.session_state.stage == 1:
        if users:
            if len(users) == 2:
                references = {"User 2": users[1]}
            st.session_state["results"] = getClosest(references, users[0], haplogroup, st.session_state.comparison_mode)
            event = show_table(table)

            best_score = st.session_state.results["Genetic Distance"].min()
            closest_names = list(st.session_state.results[st.session_state.results["Genetic Distance"] == best_score].index)
            closest_data = pd.DataFrame(get_metadata(closest_names), columns=["Name", "Country", "Group", "LAT", "LON", "Sex", "MTHG"])
            closest_data.LAT = [float(x.replace(",", ".")) for x in closest_data.LAT]
            closest_data.LON = [float(x.replace(",", ".")) for x in closest_data.LON]
        else:
            st.write("Please upload a user file or switch to 'Ancient vs. Ancient' mode.")
            set_stage(0)

        if event.selection.rows:
            st.session_state["individual_results"] = st.session_state.results.iloc[event.selection.rows[0]]
            id = st.session_state.individual_results.name

            meta = get_metadata(id)
            meta = pd.Series(meta[1:], index=metadata_header[1:], name=id)
            st.write(pd.concat([st.session_state.individual_results, meta]))

        if(st.session_state.show_map):
            st.map(closest_data)

main()
