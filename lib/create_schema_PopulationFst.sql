CREATE TABLE IF NOT EXISTS PopulationFst (
# Standard columns
label VARCHAR(100), # Dataset label
chr VARCHAR(50),
start INT, # Start site, always specified
end INT,   # End site, specified by END or assumed start + 1
name VARCHAR(100),  # Provided name
amreas_Nhat DOUBLE, # Fst Numerator: American - East Asian
amrsas_Nhat DOUBLE, # Fst Numerator: American - South Asian
amreur_Nhat DOUBLE, # Fst Numerator: American - European
amrafr_Nhat DOUBLE, # Fst Numerator: American - African
eassas_Nhat DOUBLE, # Fst Numerator: East Asian - South Asian
easeur_Nhat DOUBLE, # Fst Numerator: East Asian - European
easafr_Nhat DOUBLE, # Fst Numerator: East Asian - African
saseur_Nhat DOUBLE, # Fst Numerator: South Asian - European
sasafr_Nhat DOUBLE, # Fst Numerator: South Asian - African
eurafr_Nhat DOUBLE, # Fst Numerator: European - African
allpop_Nhat DOUBLE, # Fst Numerator: All Continental Populations
amreas_Dhat DOUBLE, # Fst Denominator: American - East Asian
amrsas_Dhat DOUBLE, # Fst Denominator: American - South Asian
amreur_Dhat DOUBLE, # Fst Denominator: American - European
amrafr_Dhat DOUBLE, # Fst Denominator: American - African
eassas_Dhat DOUBLE, # Fst Denominator: East Asian - South Asian
easeur_Dhat DOUBLE, # Fst Denominator: East Asian - European
easafr_Dhat DOUBLE, # Fst Denominator: East Asian - African
saseur_Dhat DOUBLE, # Fst Denominator: South Asian - European
sasafr_Dhat DOUBLE, # Fst Denominator: South Asian - African
eurafr_Dhat DOUBLE, # Fst Denominator: European - African
allpop_Dhat DOUBLE, # Fst Denominator: All Continental Populations
amreas_Fst DOUBLE, # Fst: American - East Asian
amrsas_Fst DOUBLE, # Fst: American - South Asian
amreur_Fst DOUBLE, # Fst: American - European
amrafr_Fst DOUBLE, # Fst: American - African
eassas_Fst DOUBLE, # Fst: East Asian - South Asian
easeur_Fst DOUBLE, # Fst: East Asian - European
easafr_Fst DOUBLE, # Fst: East Asian - African
saseur_Fst DOUBLE, # Fst: South Asian - European
sasafr_Fst DOUBLE, # Fst: South Asian - African
eurafr_Fst DOUBLE, # Fst: European - African
allpop_Fst DOUBLE, # Fst: All Continental Populations
gd_id BIGINT, # GenomicData direct-reference key
PRIMARY KEY(label,name,chr,start,end),
KEY(gd_id),
KEY(label,chr,start,end),
KEY(label,allpop_Fst)
)