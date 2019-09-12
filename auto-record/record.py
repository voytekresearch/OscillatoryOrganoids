"""
    A zmq client to test remote control of open-ephys GUI
    Documentation can be found at:
    https://open-ephys.atlassian.net/wiki/spaces/OEW/pages/23265310/Network+Events
"""

import zmq
import os
import time


def record(record_length=180):

    # Basic start/stop commands
    start_cmd = 'StartRecord'
    stop_cmd = 'StopRecord'
    # Example settings
    rec_dir = os.path.join(os.getcwd(), 'Output_RecordControl')
    print "Saving data to:", rec_dir
    # Connect network handler
    ip = '127.0.0.1'
    port = 5556
    timeout = 1.

    url = "tcp://%s:%d" % (ip, port)

    with zmq.Context() as context:
        with context.socket(zmq.REQ) as socket:
            # Set timeout in milliseconds
            socket.RCVTIMEO = int(timeout * 1000)
            # Connect to Open Ephys
            socket.connect(url)
            # Start data acquisition
            socket.send('StartAcquisition')
            answer = socket.recv()
            print answer
            time.sleep(5)
            # Begin recording of data
            socket.send(start_cmd)
            answer = socket.recv()
            print answer
            time.sleep(record_length)
            socket.send(stop_cmd)
            answer = socket.recv()
            print answer
            # Finally, stop data acquisition; it might be a good idea to 
            # wait a little bit until all data have been written to hard drive
            time.sleep(5)
            socket.send('StopAcquisition')
            answer = socket.recv()
            print answer

#run_client()
