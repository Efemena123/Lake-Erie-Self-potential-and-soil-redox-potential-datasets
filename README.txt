This package contains field data acquired at both Crane Creek (CRC) and Old Woman Creek (OWC) COMPASS sites, along with Python scripts supporting the submitted manuscript on "moisture content-driven reversals in self-potential–redox relationships in soils across the Lake Erie coastal interface zones". It includes field data collected from self-potential (SP), soil redox potential (Eh), soil electrical conductivity from the electromagnetic induction (EMI) method, in-situ soil moisture content (SMC), and soil temperature (T) datasets, including profile distance (m) and coordinates (latitude and longitude), which were acquired one day per month across land–lake interfaces along the Lake Erie Great Lake. These datasets were useful for assessing how self-potential relates to soil redox potential in the Great Lakes interface soils. We combined two geophysical methods, SP and EMI, and data from soil sensors to investigate spatial heterogeneities across land–lake interfaces.

In addition, the package contains a Spyder Python script for inverting the SP dataset. When run, the script automatically launches a GUI dashboard where you can import your CSV data and enter your grid cell (i.e., your field electrode spacing); overall, it is user-friendly.

The data archive is organized in the following directories:

1) SP_Eh_CRC

Includes two set-up profile zones with their respective distance and coordinate locations (Latitude & Longitude): wetland identified as (CCW) and transition identified as (CCT), with repeated monthly measurements (June_2023 – October_2024) provided as CSV files, alongside a text file detailing dataset units.

2) SP_Eh_OWC

Includes two set-up profile zones with their respective distance and coordinate locations (Latitude & Longitude): wetland identified as (OWC) and wetland edge identified as (OWCED), with repeated monthly measurements (June_2023 – March_2025) provided as CSV files, alongside a text file detailing dataset units.

3) SP_redox_2complete

This CSV file is set up for the Python script for inverting or reproducing the measure SP response.

4) MS2_SP_Ehcode.py

The Python script name in the file.

--> data units
SP is measured in millivolts (mV)
Eh is measured in millivolts (mV)
EC is measured in millisiemens per meter (mS/m)
SMC is expressed as a percentage (%)
T is measured in degrees Celsius (0C)


Remark: The ".csv) files can be opened with MS Excel or any other spreadsheet software.
