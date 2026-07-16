import serial
import time
import sys

# CHANGE THIS to match the port shown in your Arduino IDE
SERIAL_PORT = '/dev/cu.usbmodem101'  # Example for macOS, change as needed
BAUD_RATE = 9600

print("==============================================================")
print("VIROVORE BIOREACTOR CYBER-BIOTIC MONITORING SYSTEM")
print("==============================================================")

try:
    # Open the serial connection
    ser = serial.Serial(SERIAL_PORT, BAUD_RATE, timeout=1)
    time.sleep(2) # Allow the connection to initialize
    print(f"[*] Successfully connected to hardware on port: {SERIAL_PORT}")
    print("[*] Waiting for incoming bioreactor telemetry...\n")
    
    while True:
        if ser.in_waiting > 0:
            # Read the raw line from the USB connection
            raw_line = ser.readline().decode('utf-8').strip()
            
            # Split the data into temperature and alert status
            data_points = raw_line.split(',')
            
            if len(data_points) == 2:
                temperature = data_points[0]
                status = data_points[1]
                
                # Check for reading errors from the sensor
                if temperature == "ERROR":
                    print("[!] SENSOR ERROR: Check DHT11 wiring.")
                    continue
                
                # Display high-value visual telemetry matching your career blueprint
                if status == "ALERT":
                    # \a triggers your Mac terminal's system "bell" sound!
                    sys.stdout.write("\a") 
                    sys.stdout.flush()
                    print(f"[\033[91m🚨 ALERT\033[0m] Temp: {temperature}°C | Bioreactor Vent: \033[91mOPEN [Venting Heat]\033[0m")
                else:
                    print(f"[\033[92m🟢 SAFE\033[0m ] Temp: {temperature}°C | Bioreactor Vent: \033[92mCLOSED [Stable]\033[0m")
                    
except KeyboardInterrupt:
    print("\n[!] Telemetry feed stopped by user. Exiting safely.")
    ser.close()
    sys.exit()
except Exception as e:
    print(f"\n[X] Connection Error: {e}")
    print("[*] Double check your USB cable and make sure the Arduino Serial Monitor is closed.")