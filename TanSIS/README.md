### TanSIS Sampling Plan

1. download the data

    * create a folder called "data" in the repo folder;

    * download the zip file from

        [Dropbox TanSIS folder: TanSIS_samplingdata.zip](https://www.dropbox.com/s/ihsiaweyvc8y119/TanSIS_samplingdata.zip?dl=0)

    * unzip data.zip, and put all the data files in the "data" folder.

2. run the script:

    * pip install all the required python modules by running `pip` in the **repo directory**: 

    ```sh
        pip install -r requirement.txt
    ```

    * then run the python script:

    ```python
        python TanSIS_sampling.py
    ```
3. results are saved in output/ folder: one folder for each district, and a csv file: **sampled\_locs\_total.csv**, including all sampled locations.
