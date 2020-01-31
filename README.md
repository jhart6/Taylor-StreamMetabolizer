# Taylor-StreamMetabolizer
Stream Metabolism on Taylor Creek, Seattle WA

This repository contains resources to build a stream metabolism model for sampling sites in  Lower Taylor Creek, Seattle WA using the R package "StreamMetabolizer."  Included in this repository: StreamMetabolizer introductory scripts, time series datasets collected in situ, and R scripts to build basic metabolism models for Taylor Creek "Mouth."

Field data collection:
1. Sensor Locations
    * Sensors were located at four sites on Lower Taylor Creek, Seattle WA:     * Mouth (47.510209,	-122.24631)
      * Garden (47.51142003,	-122.2468133)
      * Hound (47.51069868,	-122.2480442)
      * Culvert (47.50948833,	-122.2481793)
2. Sensor Descriptions
    * Each location was fitted with three sensors, firmly attached to a 4ft piece of rebar.
        1. PME MiniDOT Sensor, fully submerged, equipped with copper anti-biofouling plate, sampled at 10min intervals
        2. Onset HOBO Conductivity Logger, fully submerged, sampled at 10 min intervals
        3. HOBO Pendant PAR Sensor, located just above water level, sampled at 10min intervals
3. Data Collation
    * Data were downloaded periodically over the sampling period.
    * Data were manually collated in Microsoft Excel, according to timestamp.

Each metabolism dataset contains the following variables:
1. Epoch.Sec: This is the unit in which MiniDOTs record DateTime.
2. DateTime.sensor: Epoch.Sec was converted to DateTime.sensor using the below equation in Excel.
  * =(Epoch.Sec/(24*60*60) + DATE(1970,1,1))
3. DateTime.adj: Epoch.Sec on the MiniDOTs is not accurate; it records a DateTime stamp that is approximately three hours ahead. I have been correcting this using the below equation in Excel.
  * = DateTime.sensor - 0.3
4. BV.volts
5. Temp.C: Recorded by PME MiniDOT concurrently with DO.
6. DO.mgL: Recorded by PME MiniDOT.
7. Q
8. Intensity.lumft2: Recorded by HOBO Pendant PAR sensor.

Dataset Details:
* Sensors were initially deployed at all four locations on June 6, 2019.
* Data were downloaded periodically over the next 6-7 months.
* Final data download performed by J. Hart on January 3, 2020.
* No reliable miniDOT data for Garden, Hound, & Culvert from January 3, 2020 to February 7, 2020.
* There is a period of poor oxygen data at the mouth site in July 2019. These data points have NOT yet been removed from the time series, but the metabolism model would benefit from their removal. Consensus is that these numbers were a fluke.
* Water Temperature time series was recorded concurrently with Oxygen on the MiniDOTs.
* According to bucket experiments in the lab, each miniDOT consistently underestimates DO concentrations (mg/L and % saturation), according to the YSI EXO2, which was calibrated for DO before every field deployment. A correction factor has been calculated for each miniDOT and should be applied to the metabolism datasets before proceeding with metabolism model.
  * Mouth
    * mg/L offset:
    * % saturation offset:
  * Garden
    * mg/L offset:
    * % saturation offset:
  * Hound
    * mg/L offset:
    * % saturation offset:
  * Culvert
    * mg/L offset:
    * % saturation offset:

Suggested Next Steps:
* Correct miniDOT time series data according to documented correction factor.
* Starting with the mouth site, choose a model setup and calibrate it well for 3 days. 
* Apply that model calibration to the entire time series.
* Apply model calibration (if appropriate) to the other three sites on Lower Taylor Creek.
