
import csv
import paramiko
import os
from datetime import datetime
import time as Time

class Logger:
    def __init__(self):
        self.filename = None
        self.startTime = Time.time()
        self.bmi0TimeOffset = None
        self.bmi1TimeOffset = None
        
    def uptime(self):
        return Time.time() - self.startTime
    
    def _connect_ssh(self):
        self.__ssh = paramiko.SSHClient()
        self.__ssh.set_missing_host_key_policy(paramiko.AutoAddPolicy())
        self.__ssh.connect(
            "192.168.3.4",
            username=os.getenv("SSH_USERNAME"),
            password=os.getenv("SSH_PASSWORD"),
        )
        print("Connected to SSH")
        
    def _upload_ssh(self, filename):
        self._connect_ssh()
        sftp = self.__ssh.open_sftp()
        sftp.put(
            filename,
            f"/home/{os.getenv('SSH_USERNAME')}/Documents/MATLAB/{filename.split('/')[-1]}",
        )
        sftp.close()
        self.__ssh.close()
            
    def start_log(self):
        self.filename = f"./measurement_log_{str(datetime.now().strftime('%Y-%m-%d %H:%M:%S'))}.csv"
        self.csvfile = open(self.filename, "w", newline="")
        # Create a CSV writer object
        self.writer = csv.writer(self.csvfile)
        self.log_row(["TIME", "A_X_0", "A_Y_0", "A_Z_0", "A_X_1", "A_Y_1", "A_Z_1", "VALVE"])
    
    def log_row(self, row):
        self.writer.writerow(row)
    
    def bmi_time_sync(self, deviceid, bmiTime):
        if deviceid==0:
            self.bmi0TimeOffset =  bmiTime - self.uptime()
        elif deviceid==1:
            self.bmi1TimeOffset = bmiTime - self.uptime()
        else:
            raise ValueError("Invalid device id")
    
    def log_bmi(self, deviceid, read_dict):
        if deviceid==0:
            logRow = [self.uptime(), read_dict['A_X'], read_dict['A_Y'], read_dict['A_Z'] , None, None, None, None]
        elif deviceid==1:
            logRow = [self.uptime(), None, None, None,  read_dict['A_X'], read_dict['A_Y'], read_dict['A_Z'], None]
        else:
            raise ValueError("Invalid device id")
        self.log_row(logRow)
    
    def log_valve(self, position):
        self.log_row([self.uptime(), None, None, None, None, None, None, position])
    
    def end_log(self, upload=True, delete=False):
        self.csvfile.close()
        if upload:
            self._upload_ssh(self.filename)
    