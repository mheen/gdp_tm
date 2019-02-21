# gdp_tm
Deriving transport matrices from observed drifter locations from NOAA's Global Drifter Program

## Getting started
Head to http://www.aoml.noaa.gov/envids/gld/FtpInterpolatedInstructions.php and download NOAA's Global Drifter Program (GDP) drifter data.
You'll need both the data (buoydata_xxx.dat) and the meta data (dirfl_xxx.dat) files.

In the gdp_tm/drifters folder create a file "dirs.json". This file points to relevant file locations that "read_gdp_data" uses.
It should contain the following fields.
```
{
 "dirs":
  {
   "input_dir": "my_directory/with_gdp/original_data/",
   "output_path": "my_directory/with_gdp/matlab_data.mat",
   "data_filenames": ["buoydata_1_5000.dat","buoydata_5001_10000.dat"],
   "meta_data_filenames": ["dirfl_1_5000.dat","dirfl_5001_10000.dat"]
  }
 }
 ```
 
 Now run "gdp_data = read_gdp_data()". This will load all the GDP drifter data into a Matlab data structure
 (and save it to the .mat file you specified in "dirs.json" so you don't have to do this again all the time).
 
