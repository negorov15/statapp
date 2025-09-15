import pandas as pd
from ete3 import NCBITaxa

ncbi = NCBITaxa()
#ncbi.update_taxonomy_database()

def get_lineage_from_cell(taxa_name):

    # Step 1: Convert name to taxid
    name2taxid = ncbi.get_name_translator([taxa_name])
    taxid = name2taxid[taxa_name][0]

    # Step 2: Get full lineage taxids
    lineage = ncbi.get_lineage(taxid)

    # Step 3: Map lineage taxids to names
    names = ncbi.get_taxid_translator(lineage)

    # Step 4: (Optional) Get ranks for each taxid
    template = ["Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"]
    ranks = ncbi.get_rank(lineage)

    # Build ordered path (Domain → ... → species/order/etc.)
    taxonomy_path = [names[t] for t in lineage if t in template]

    return taxonomy_path

def get_lineage(dataset):
    """
        Reads a tab-delimited dataset containing a column 'Taxa' with NCBI taxonomy IDs
        and returns their corresponding full taxonomic lineages.

        Steps:
        1. Reads the dataset into a DataFrame.
        2. Iterates through the Taxa IDs:
            - If the ID is -2, appends 'Not assigned' to the result list.
            - Otherwise, queries the ETE3 NCBITaxa database to fetch the full lineage.
            - Skips invalid/missing taxon IDs that raise ValueError.
        3. Passes the collected lineages through:
            - sort_ranks(): to arrange lineage elements consistently by rank.
            - translate_lineage(): to convert taxon IDs into human-readable taxonomic names.
        4. Returns a list of translated lineage strings.

        Parameters
        ----------
        dataset : str
            Path to a tab-delimited file (.tsv) with at least one column named 'Taxa'.

        Returns
        -------
        list
            A list of strings where each string represents the full lineage of a taxon.
            Example: ["Bacteria; Proteobacteria; Gammaproteobacteria; Enterobacterales; ..."]
            For -2 entries, "Not assigned" is returned instead.
        """

    # Read the input file
    df = pd.read_csv(dataset, delimiter="\t")
    list_of_lineages = []
    # Extract taxa ids and loop through the list
    tax_ids = df['Taxa']
    for id in tax_ids:
        # -2 = 'Not assigned'.
        # We just add 'not assigned' rank also to the list to avoid misunderstandings
        if id == -2:
            list_of_lineages.append('Not assigned')
            continue
        else:
            try:
                list_of_lineages.append(ncbi.get_lineage(id))
            except ValueError:
                continue
    # Call 'sort_ranks' function and store it to the variable 'sorted_lineages'
    sorted_lineages = sort_ranks(list_of_lineages)
    return translate_lineage(sorted_lineages)

# Helper function. Sorts ranks so that only major ranks left.
# Input: List of lineages
# Output: List of sorted lineages
def sort_ranks(data):
    major_ranks = ["superkingdom", "domain", "phylum", "class", "order", "family", "genus", "species"]
    sorted_ranks = []
    for lineage in data:
        # Temporary list for storing ranks.
        # Is required, so the final list would have separate lists of ranks associated with each element
        lineage_ranks = []
        # Check if the element of list has an actual lineage.
        # If it is the case, just add 'Not assigned' to the list
        if lineage == "Not assigned":
            sorted_ranks.append([lineage])
        else:
            # Get ranks for each lineage
            rank = ncbi.get_rank(lineage)
            for i in range(len(lineage)):
                # If the rank is not in the list of major ranks, we skip it
                if rank[lineage[i]] in major_ranks:
                    lineage_ranks.append(lineage[i])
                else:
                    continue
            # add the temporary list to the main list
            sorted_ranks.append((lineage_ranks))
    return sorted_ranks

# Function translates the lineages into names.
# Input: List of lineages
# Output: Translated list of lineages
def translate_lineage(list_of_lineages):
    # Empty list where all names will be stored
    names = []
    # Loop through each lineage list and get the names of each rank
    for lineage in list_of_lineages:
        # If the rank is 'not assigned' simply add the value to the list.
        if lineage == ["Not assigned"]:
            names.append(lineage)
        else:
            # get_taxid_translator(lineage) creates a dictionary,
            # so we need to get all values and store it to the list
            translate = ncbi.get_taxid_translator(lineage)
            names.append([translate[taxid] for taxid in lineage])
    # Create a DataFrame where:
    # Rows: taxa
    # Columns: taxonomy lineage associated with the taxa
    df = pd.DataFrame(names, columns=["Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"])
    return df