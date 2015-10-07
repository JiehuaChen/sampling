##Ghana Sampling Plan


1. clone the data 
2. download the data

	* create a folder called "data" in the repo folder;
	* download both geosurvey prediction data and Ghana shapefile from
	
	[geosurvey data](https://www.dropbox.com/s/f0fyy3dlh44ar0b/geosurvey_results.zip?dl=0)
	
	[shapefile](https://www.dropbox.com/s/fj61a0rvxyji23d/shapefile.zip?dl=0)


	* unzip the zip files, and put all the data files in the "data" folder.

3. run the script:
	
	pip install all the required python modules by running pip in the repo directory:
    `pip install -r requirement.txt`
then run the python script:
    `python TanSIS_sampling.py`

4. results are saved in output/ folder: one folder for each district, and a csv file: sampled_locs_total.csv, including all sampled locations.
