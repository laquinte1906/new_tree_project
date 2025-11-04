README.md file was generated on 2024-05-19 by Ian Quintas


GENERAL INFORMATION

1. Title of Dataset: Data from: Seasonal effects of farmer-managed livestock grazing exclusions on bird communities in Burkina Faso

2. Author Information
	Corresponding Author 
		Name: Dr Gabriel Marcacci
                Institution: Swiss Ornithological Institute, Sempach, Switzerland
                Email: gabriel.marcacci@vogelwarte.ch

	Co-author 1
		Name: Mr Ian Quintas
		Institution1: Swiss Ornithological Institute, Sempach, Switzerland
		Institution2: Department of Ecology and Evolution, University of Lausanne, Switzerland
        
	Co-author 2
		Name: Mr Ambroise N Zongo
		Insitution: tiipaalga, Ouagadougou, Burkina Faso

	Co-author 3
		Name: Dr Pius Korner
                Institution: Functional Swiss Ornithological Institute, Sempach, Switzerland
	
        Co-author 4
		Name: Mrs Alexandra Kuttnig
		Institution1: Swiss Ornithological Institute, Sempach, Switzerland
		Institution2: Department of Environmental Sciences, University of Basel, Switzerland

	Co-author 5
		Name: Dr Reto Spaar
		Institution: Swiss Ornithological Institute, Sempach, Switzerland

	Co-author 6
		Name: Mr Bakary Diakité
		Insitution: tiipaalga, Ouagadougou, Burkina Faso

        Co-author 7
		Name: Mrs Franziska Kaguembèga-Müller
		Institution1: tiipaalga. Ouagadougou, Burkina Faso
		Institution2: newTree, Bern, Switzerland
		
        Co-author 8
		Name: Dr Alain Jacot
		Insitution: Swiss Ornithological Institute, Sempach, Switzerland
		




3. Date of data collection: November 2022 - March 2023

4. Geographic location of data collection: Central Burkina Faso

5. Funding sources: NA

6. Recommended citation for this dataset: Quintas, Marcacci et al. (2025), Data from: Seasonal effects of farmer-managed livestock grazing exclusions on bird communities in Burkina Faso 



DATA & FILE OVERVIEW

1. Data collection

Birds: Birds were surveyed via Passive Acoustic Monitoring with Autonomous Recording Units (Audiomoths) in 27 landscape units, each including three study sites (exclosure, open control, woody control) in 2 survey rounds (dry and wet). 
Initially N = 162 bird surveys but four sites were removed because of technical issues with the acoustic recorders (N = 158 bird surveys).

Vegetation: 
In each site all tree species and their number of stems were recorded within a radius of 30 meters around the ARU for four vegetation strata (below 1 meter, between 1 and 2 meters, between 2 and 5 meters and above 5 meters).
The coverage and mean height of herbaceous vegetation was estimated within a 10-meters radius around the ARU. 
The tree inventory was only conducted once, whereas the herbaceous vegetation variables were measured for each season to account for temporal changes.

NDVI: 
Remote sensing was used to calculate vegetation cover inside every exclosure and control site (within a 50-meters radius) but also for each 500-meters radius landscape unit, informing about the surrounding landscape vegetation.
The NDVI was calculated from 5 satellite images (Sentinel-2 with a 10-meters resolution) for each season and the mean, median and maximum values were computed.

2. File List: 
Notes: 
"grazing exclusion" (or "exclosure") are denoted as "enclos" in the dataframes. 
"wet" season is denoted as "post rainy" season in the dataframes.

	 1.   Name: by_survey.csv
	      Description: Dataset containing all variables summarized by study sites and seasons (i.e., surveys, N = 158).
	    
	 3.   Name: t.dat0.csv
	      Description: All detections of species (number of 1-min recordings each species was detected in a given survey) per survey. Only species with at least 30 detections in total were considered for this analysis.
	      
	 4.   Name: all_species.csv
	      Description: Functional and life-history traits for every species recorded in the study.
	      

DATA SPECIFIC INFORMATION

1. by_survey.csv

landscape_id: landscape unit id (factor)
year: year in which the grazing exclusion was established (factor)
age_enclos: age of each grazing exclusion at the time of data collection (continuous)
habitat: exclosure, open or woody (factor with 3 levels)
site_id: name of the study sites (=landscape_id + habitat) (factor)
season: type of season (dry or post rainy (= wet)) (factor with 2 levels)
season_site_id: identifyer for season + site_id (factor)
bird_rich: number of bird species (continuous)
dec_long: longitude (in decimal degree)
dec_lat: latitude (in decimal degree)
tot_strata: Total number of trees across all strata (continuous)
plant_rich_tot: Total number of plant species across all strata (continuous)	
Shannon_vertical: vertical vegetation heterogeneity (Shannon index across strata) (continuous) 	
herb_cov: coverage of herbaceous vegetation within a 10-meters radius (continuous)
herb_height: mean height of herbaceous vegetation within a 10-meters radius (in meters) (continuous)
max.NDVI.50: max NDVI within a 50-meters radius	(site-level) (continuous)
max.NDVI.500: max NDVI within a 500-meters radius (landscape-level) (continuous)	
nb_houses: number of houses per landscape unit used as proxy for grazing pressure (continuous)	
diff.max.NDVI.500: Difference of max NDVI between seasons for each landscape unit (continuous)


2. t.dat0.csv
species_sci: scientific species names (character)
habitat: exclosure, open or woody (factor with 3 levels)
Npres: number of detections of each species (=number of 1-min recordings in which each species were detected) (from 1-20; continuous)
Nabs: number of absences of each species (=number of 1-min recordings in which each species were not detected) (from 1-20; continuous)
landscape_id: id of the landscape unit (factor)
season: type of season (dry or post rainy (= wet)) (factor with 2 levels)
Nabs.interv: number of absences within species-specific time intervals (interval from the first to the last 1-minute recording having at least a minimal probability for detection of >0.001) (from 1-20; continuous)


3. all_species.csv
species_sci: scientific species names (character)
habitat_pref: preferred habitat of the species (open/intermediate/forest) (factor with 3 levels)
migration: migration strategy adopted by the species (resident/migrant) (factor with 2 levels)
diet: preferred diet of the species (factor with 6 levels)









