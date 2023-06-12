from Bio import Entrez, Medline
from main import *

refs = rec.annotations['references']

for ref in refs:
  if ref.pubmed_id != '':
    print(ref.pubmed_id)

handle = Entrez.efetch(db='pubmed', id=[ref.pubmed_id], rettype='medline', retmode='text')
records = Medline.parse(handle)

for med_rec in records:
   print(med_rec)