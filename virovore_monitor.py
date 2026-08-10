import serial
import time
import sys
import math
import csv
import os
from contact_scorer import calculate_hydrophobic_score

# CHANGE THIS to match the port shown in your Arduino IDE
SERIAL_PORT = '/dev/cu.usbmodem101' 
BAUD_RATE = 9600
PDB_PATH = "docking.res.1.pdb"
LOG_FILE = "telemetry_log.csv"

# Initialize/Clear the CSV Log File on startup
with open(LOG_FILE, mode='w', newline='') as f:
    writer = csv.writer(f)
    # Write the header columns
    writer.writerow(["Timestamp", "Temperature", "HydrophobicScore", "Contacts", "Status"])

print("==============================================================")
print("VIROVORE BIOREACTOR LOGGING SYSTEM ACTIVED")
print("==============================================================")

# Read the baseline structural data
try:
    print("[*] Launching Silicon Virovore pipeline parser...")
    baseline_score, baseline_contacts = calculate_hydrophobic_score(PDB_PATH)
    print(f"[+] Pipeline verified baseline structural data:")
    print(f"    - Baseline Hydrophobic Contacts: {baseline_contacts}")
    print(f"    - Baseline Hydrophobic Score: {baseline_score}\n")
except Exception as e:
    print(f"[X] Pipeline Error: Could not read {PDB_PATH}. Using fallback defaults. Error: {e}")
    baseline_score, baseline_contacts = 175.03, 222

try:
    ser = serial.Serial(SERIAL_PORT, BAUD_RATE, timeout=1)
    time.sleep(2) # Allow connection to settle
    print(f"[*] Connected to Elegoo hardware on: {SERIAL_PORT}")
    print("[*] Logging and streaming telemetry to CSV...\n")
    
    start_time = time.time()
    
    while True:
        if ser.in_waiting > 0:
            raw_line = ser.readline().decode('utf-8').strip()
            data_points = raw_line.split(',')
            
            if len(data_points) == 2:
                try:
                    temp = float(data_points[0])
                except ValueError:
                    continue
                status = data_points[1]
                
                # Calculate dynamic thermal denaturation
                if temp > 25.0:
                    decay_factor = math.exp(-0.15 * (temp - 25.0))
                    current_score = baseline_score * decay_factor
                    current_contacts = int(baseline_contacts * decay_factor)
                else:
                    current_score = baseline_score
                    current_contacts = baseline_contacts
                
                elapsed_time = round(time.time() - start_time, 1)
                
                # Append data to CSV Log file for Streamlit to read
                with open(LOG_FILE, mode='a', newline='') as f:
                    writer = csv.writer(f)
                    writer.writerow([elapsed_time, temp, round(current_score, 2), current_contacts, status])
                
                # Console display
                if status == "ALERT":
                    sys.stdout.write("\a") 
                    sys.stdout.flush()
                    print(f"[🚨 ALERT] {elapsed_time}s | Temp: {temp:.2f}°C | Score: {current_score:.2f} | Contacts: {current_contacts}")
                else:
                    print(f"[🟢 SAFE ] {elapsed_time}s | Temp: {temp:.2f}°C | Score: {current_score:.2f} | Contacts: {current_contacts}")
                          
except KeyboardInterrupt:
    print("\n[!] Logging stopped. CSV saved safely.")
    ser.close()
    sys.exit()
except Exception as e:
    print(f"\n[X] Connection Error: {e}")