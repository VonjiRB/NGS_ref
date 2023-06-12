# Entrez: provides access to the NCBI's Entrez utilities
# SeqIO: allows parsing and manipulating sequence
from Bio import Entrez, SeqIO, Medline
Entrez.email = 'vonjirabe40@gmail.com'

# to search the nucleotide database for records matching the given criteria
# the gene name is "CRT" and the organism is "Plasmodium falciparum"
handle = Entrez.esearch(db='nucleotide', term='CRT [Gene Name] AND "Plasmodium falciparum"[Organism]')
# The search results are read using Entrez.read and stored in rec_list. 
rec_list = Entrez.read(handle)

# These lines check if the number of retrieved records (rec_list['RetMax']) 
# is less than the total number of records matching the search criteria (rec_list['Count'])
# If there are more records than initially retrieved, it performs a new search 
# with retmax parameter set to the total count to retrieve all the records.
if rec_list['RetMax'] < rec_list['Count']:
    handle = Entrez.esearch(db='nucleotide', term='CRT [Gene Name] AND "Plasmodium falciparum"[Organism]',
                            retmax = rec_list['Count'])
rec_list = Entrez.read(handle)

# This line extracts the list of record identifiers (IDs) from the search results, which will be used
# to retrieve the corresponding sequences.
id_list = rec_list['IdList']

# These lines use the Entrez.efetch function to retrieve the nucleotide sequences corresponding to 
# the given list of IDs.
# The sequences are fetched in GenBank format (rettype='gb'). The SeqIO.parse function from SeqIO 
# module is used to parse the fetched sequences and store them in the recs list.
hdl = Entrez.efetch(db='nucleotide', id=id_list, rettype='gb')
recs = list(SeqIO.parse(hdl, 'gb'))

# This loop iterates over the fetched records in recs. It checks if the record's name matches 'KM288867'. 
# If found, it breaks the loop and proceeds with the subsequent code.
for rec in recs:
    if rec.name == 'KM288867':
        break

# This loop iterates over the features of the rec record. Features are specific regions or elements within
# a sequence, such as genes, exons, promoters, etc.
for feature in rec.features:
    # If the feature type is identified as a 'gene', the code retrieves the value of the 'gene' qualifier 
    # from the feature and prints it.
    # The 'gene' qualifier typically provides the name or identifier of the gene associated with the feature.
    if feature.type == 'gene':
        print(feature.qualifiers['gene'])

    # If the feature type is identified as an 'exon', the code retrieves the location information of the
    # feature (start position, end position, and strand) and prints them.
    # The location information is stored in the location attribute of the feature.
    
    # If the feature type is 'gene', it prints the value of the 'gene' qualifier associated with that feature
    # (feature.qualifiers['gene']).

    # If the feature type is 'exon', it retrieves the location information of the feature (feature.location) and
    # prints the start position (loc.start), end position (loc.end), and strand information (loc.strand).

    # If the feature type is neither 'gene' nor 'exon', it simply prints a message indicating that the feature 
    # was not processed.

    elif feature.type == 'exon':
        loc = feature.location
        print(loc.start, loc.end, loc.strand)
    else:
        print('not processed:\n%s' % feature)

# This loop iterates over the annotations of the record (rec.annotations). It prints each annotation's name and value 
# using a formatted string.
              
for name, value in rec.annotations.items():
    print('%s=%s' % (name, value))

# This line extracts the sequence from the record (rec.seq) and assigns it to the sequence variable for further 
# processing or analysis.
sequence = rec.seq

# In summary, these lines of code process the features and annotations of a specific record obtained from the 
# fetched records. It prints information about genes and exons, as well as any other features that are not of 
# the gene or exon type. It also prints the annotations associated with the record and extracts the sequence 
# for further use.






