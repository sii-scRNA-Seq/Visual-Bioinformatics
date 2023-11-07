from flask import Flask, jsonify
from flask_cors import CORS

app = Flask(__name__)
CORS(app)

@app.route('/')
def get_data():
    data = {
        'title': 'Hello, World!',
        'message': "You can't see this bit"
    }
    return jsonify(data)

if __name__ == '__main__':
    app.run()