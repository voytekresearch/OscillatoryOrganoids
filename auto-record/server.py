from flask import Flask, render_template, Markup, Response
from flask_socketio import SocketIO, emit
from multiprocessing import Process
from record import run_client
import time
import os
import sys
import argparse

url = 'localhost:5000'

app = Flask(__name__)
app.config['SEND_FILE_MAX_AGE_DEFAULT'] = 0
app.config['SECRET_KEY'] = 'secret!'
socketio = SocketIO(app)

myFile=args.file
fileLocation=myFile

@app.route('/')
def home():
    response=Response(render_template('main.html', jsloc=jsloc, css=css, body=Markup(body), stderr=Markup(stderr), stdout=Markup(stdout))
)
    response.headers["Cache-Control"] = "no-cache, no-store, must-revalidate" # HTTP 1.1.
    response.headers["Pragma"] = "no-cache" # HTTP 1.0.
    response.headers["Expires"] = "0" # Proxies
    return response
@socketio.on('connect')
def connect():
    print('Connected to client')
@socketio.on('disconnect')
def disconnect():
    print('Client disconnected')
def main():
    schedule.every().day.at('14:30').do(run_client)
    def startFlask():
        socketio.run(app, port=5000, debug=True)
    try:
        flask=Process(target=startFlask)
        flask.start()
    except KeyboardInterrupt:
        print("Terminating flask server.")
        flask.terminate()
        flask.join()
    while True:
        schedule.run_pending()
        time.sleep(1)
