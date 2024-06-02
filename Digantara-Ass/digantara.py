from datetime import datetime, timedelta
from sgp4.api import Satrec, jday
from skyfield.api import EarthSatellite, load
import math
from sgp4.earth_gravity import wgs84
import numpy as np
import warnings

### i am referrirng object as satellite

class Satellite:

    #this object creates a satellite class which can provide it e,r,v vectors and check the sunlit condition

    def __init__(self, line1, line2):
        self.satellite = Satrec.twoline2rv(line1, line2)

    def is_sunlit(self, curr_time): #function for ckecing sunlit condition using skyfield.api
        eph = load('de421.bsp')
        ts = load.timescale()
        t = ts.utc(curr_time.year, curr_time.month, curr_time.day, curr_time.hour, curr_time.minute, curr_time.second)
        sat = EarthSatellite(line1,line2)
        self.sunlit = sat.at(t).is_sunlit(eph)
        return self.sunlit

    def propagate(self, jd, fr):
        e, r, v = self.satellite.sgp4(jd, fr)
        return e, r, v

class Tracker:
    #tracker object that can provide e,r,v of tracker, given the TLE of trakcer

    def __init__(self, line1, line2):
        self.tracker = Satrec.twoline2rv(line1, line2)

    def propagate(self, jd, fr):
        e, r, v = self.tracker.sgp4(jd, fr)
        return e, r, v

class DetectionLogic:
    #class for forming the detection logic taking input of r,v vectors of tracker and satellite
    #i.e. when angle between the vel vector of tracker and track2sat_vec is < 30/2 deg, when the norm of track2sat_vec < 1000 Km

    def is_crossing(r_sat, v_sat, r_track, v_track, view_angle, min_dist):

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
        
class TrackerData:
    #class to create the TLE of the tracker through the given values

    def __init__(self, **kwargs):
        self.utc_epoch = utc_epoch
        self.SMA = SMA
        self.inclination = inclination
        self.eccentricity = eccentricity
        self.RAAN = RAAN
        self.arg_of_perigee = arg_of_perigee
        self.mean_anomaly = mean_anomaly
        self.coord_system = coord_system

    def calculate_checksum(self, line): #calculating checksum value to add at last of line
        checksum = 0
        for char in line:
            if char.isdigit():
                checksum += int(char)
            elif char == '-':
                checksum += 1
        return str(checksum % 10)

    def get_tracker_TLE(self):
        epoch_str = self.utc_epoch.strftime("%y%j")
        seconds_since_midnight = (self.utc_epoch - self.utc_epoch.replace(hour=0, minute=0, second=0, microsecond=0)).total_seconds()
        fraction_of_day = seconds_since_midnight / (24 * 60 * 60)
        fraction_of_day_str = f"{fraction_of_day:.8f}".lstrip("0")

        # Find time period in rev/s
        #T = 2 * math.pi * math.sqrt(self.SMA  ** 3 / (wgs84.mu))
        n = math.sqrt(wgs84.mu / self.SMA ** 3) *24*60*60/ (2 * math.pi)

        line1 = f"1 00000U 00000A   {epoch_str}{fraction_of_day_str}  .00000000  00000-0  00000-0 0  999"
        line2 = f"2 00000 {self.inclination:8.4f} {self.RAAN:8.4f} {self.eccentricity:07.0f} {self.arg_of_perigee:8.4f} {self.mean_anomaly:08.4f} {n:11.8f} 9999"

        line1 = line1 + self.calculate_checksum(line1)
        line2 = line2 + self.calculate_checksum(line2)

        return line1, line2

class Main:
    #main class to integrate all classes

    def __init__(self, satellite, tracker):
        self.satellite = satellite
        self.tracker = tracker
        

    def run_detection(self, start_time, time_step, duration_days,view_angle,min_dist2track): 
        endtime = start_time + timedelta(days=duration_days)
        detection_flag = False

        while start_time < endtime:
            if start_time.minute == 0 and start_time.second == 0:
                print("now running for hour ",start_time	,'\n')

            jd, fr = jday(start_time.year, start_time.month, start_time.day, 
                          start_time.hour, start_time.minute, start_time.second)
            
            e_Sat, r_sat, v_sat = self.satellite.propagate(jd, fr)
            e_track, r_track, v_track = self.tracker.propagate(jd, fr)

            if e_Sat > 0.001 or e_track > 0.001:
                warnings.warn(('e_Sat:{},e_track:{} error is too large').format(e_Sat,e_track))
                start_time += timedelta(seconds=time_step)
                continue

            if self.satellite.is_sunlit(start_time):
                if DetectionLogic.is_crossing(r_sat, v_sat, r_track, v_track, view_angle, min_dist2track):
                    if not detection_flag:
                        print("Eureka!!!  The object detected, signal started at ", start_time)
                    detection_flag = True
                else:
                    if detection_flag:
                        print(":(:(:(:(:(:( Signal lost at ", start_time,'\n')
                    detection_flag = False

            start_time += timedelta(seconds=time_step)

if __name__ == "__main__":
    with open('TLE.txt[2119].txt', 'r') as file:
        lines = file.readlines()
    
    line1 = lines[0]
    line2 = lines[1]

    utc_epoch = datetime(2023, 3, 21, 0, 0, 0)  # Epoch time in UTC
    # ist_epoch = utc_epoch.astimezone(pytz.timezone('Asia/Kolkata'))
    SMA = 6878  
    inclination = 97.4  
    eccentricity = 0  
    RAAN = 269.8035  
    arg_of_perigee = 331.7425  
    mean_anomaly = 0  
    coord_system = "TEME"  

    tracker_data = TrackerData(utc_epoch=utc_epoch, SMA=SMA, inclination=inclination, eccentricity=eccentricity,
                               RAAN=RAAN, arg_of_perigee=arg_of_perigee, mean_anomaly=mean_anomaly, coord_system=coord_system)

    satellite = Satellite(line1, line2)
    tracker_line1, tracker_line2 = tracker_data.get_tracker_TLE()
    tracker = Tracker(tracker_line1, tracker_line2)

    start_time = datetime(2023, 3, 21, 0, 0, 0)
    time_step = 1 # in sec
    duration_days = 1
    min_dist2track = 1000 #km
    view_angle = 30 #in deg

    main = Main(satellite, tracker)
    main.run_detection(start_time, time_step, duration_days,view_angle=view_angle,min_dist2track=min_dist2track)
