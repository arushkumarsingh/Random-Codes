from datetime import datetime, timedelta
from sgp4.api import Satrec, jday
from skyfield.api import load, Topos, utc, EarthSatellite
from sgp4.io import twoline2rv
from sgp4.earth_gravity import wgs84
import numpy as np
import TrackerData2TLE
import sys
import warnings
import math

def sunlit(line1,line2,curr_time):
    # Load the ephemeris data
    eph = load('de421.bsp')
    ts = load.timescale()
    # curr_time = datetime.utcnow()
    t = ts.utc(curr_time.year, curr_time.month, curr_time.day, curr_time.hour, curr_time.minute, curr_time.second)
    satellite = EarthSatellite(line1,line2)
    sunny = satellite.at(t).is_sunlit(eph)

    return sunny

def is_crossing(r_sat, v_sat, r_track, v_track,view_angle,min_dist):

    r_sat = np.array(r_sat)
    v_sat = np.array(v_sat)
    r_track  = np.array(r_track)
    v_track = np.array(v_track)

    track2sat_Vec = r_sat - r_track
    norm_v_sat = np.linalg.norm(v_sat)
    norm_track2sat_Vec = np.linalg.norm(track2sat_Vec)
    # print(norm_track2sat_Vec)

    angle = np.degrees(np.arccos(np.dot(v_sat, track2sat_Vec) / (norm_v_sat * norm_track2sat_Vec)))
    
    if angle < view_angle/2 and norm_track2sat_Vec < min_dist:
        return True
    else:
        return False

    

with open('TLE.txt[2119].txt', 'r') as file:
    lines = file.readlines()

line1 = lines[0]
line2 = lines[1]
tracker_line1,tracker_line2 = TrackerData2TLE.get_tracker_TLE()



satellite = Satrec.twoline2rv(line1, line2)
traker = Satrec.twoline2rv(tracker_line1,tracker_line2)

start_time = datetime(2023, 3, 21, 0, 0, 0)
time_step = 10 #in sec
duration_days = 1  
endtime = start_time + timedelta(days=duration_days)


detection_flag = 0
while start_time < endtime:
    # print("running for ",start_time,'\n')

    jd, fr  = jday(start_time.year, start_time.month, start_time.day, start_time.hour, start_time.minute, start_time.second)
    e_Sat, r_sat, v_sat = satellite.sgp4(jd, fr) #error e, position vector r and velocity vector v all in Cartesian coordinates.
    e_track, r_track, v_track = traker.sgp4(jd, fr)

    if e_Sat > 0.001 or e_track>0.001:
        warnings.warn(('e_Sat:{},e_track:{} error is too large').format(e_Sat,e_track))
        start_time += timedelta(seconds=time_step)
        continue

    if sunlit(line1,line2,start_time):
        if is_crossing( r_sat, v_sat, r_track, v_track,view_angle=30,min_dist=1000):
            if detection_flag == 0: 
                print("The object detected, signal started at ", start_time)
            detection_flag = 1
        else:
            if detection_flag == 1:
                print("Signal lost at ", start_time)
            detection_flag = 0


    start_time += timedelta(seconds=time_step)

   


    # print('TLE Satellite:',line1,line2,'\n')
    # print('TLE Tracker:',tracker_line1,'\n',tracker_line2,'\n')
    # print(('jd:{},fr: {}').format(jd,fr))
    # print(('satelllite e: {}, r:{}, v:{}').format(e_Sat,r_sat,v_sat))
    # print(('Tracker e: {}, r:{}, v:{}').format(e_track,r_track,v_track))
