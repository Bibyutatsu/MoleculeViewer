import csv
from flask import Flask, render_template, request, redirect, url_for
import requests
from pager import Pager
from rdkit import Chem
from rdkit.Chem import Draw, Descriptors
import pandas as pd
import os

APPNAME = "Molecule Viewer"
STATIC_FOLDER = 'example'
IMAGE_FOLDER = 'example/images'
TABLE_FILE = "example/fakecatalog.csv"

pd.DataFrame([], columns=['name', 'ra', 'dec', 'notes']).to_csv(TABLE_FILE, index=False)


def read_table(url):
    """Return a list of dict"""
    # r = requests.get(url)
    with open(url) as f:
        return [row for row in csv.DictReader(f.readlines())]

def add_table(url, smiles):
    df = pd.read_csv(url)
    num_of_mols = len(df)
    
    mol = Chem.MolFromSmiles(smiles)
    ra = Descriptors.ExactMolWt(mol)
    dec = Descriptors.ExactMolWt(mol)
    notes = 'Mol' + str(num_of_mols + 1)
    name = 'Mol' + str(num_of_mols + 1) + '.jpg'
    Draw.MolToImageFile(mol, os.path.join(IMAGE_FOLDER, name))
    row_to_add = pd.DataFrame([[notes, ra, dec, notes]], columns=['name', 'ra', 'dec', 'notes'])
    df = df.append(row_to_add)
    df.to_csv(url, index=False)


table = read_table(TABLE_FILE)
pager = Pager(len(table))


app = Flask(__name__, static_folder=STATIC_FOLDER)
app.config.update(
    APPNAME=APPNAME,
)


@app.route('/')
def index():
    global pager, table
    if pager.count == 0:
        add_table(TABLE_FILE, 'CC')
        table = read_table(TABLE_FILE)
        pager = Pager(len(table))
    return redirect('/0')


@app.route('/<int:ind>/')
def image_view(ind=None):
    global table, pager
    if ind >= pager.count:
        return render_template("404.html"), 404
    else:
        pager.current = ind
        return render_template(
            'imageview.html',
            index=ind,
            pager=pager,
            data=table[ind])


@app.route('/goto', methods=['POST', 'GET'])
def goto():
    return redirect('/' + request.form['index'])


@app.route('/molweight', methods=['POST', 'GET'])
def molweight():
    global table, pager
    add_table(TABLE_FILE, request.form['smiles'])
    table = read_table(TABLE_FILE)
    pager = Pager(len(table))
    return redirect('/0')


if __name__ == '__main__':
    app.run(debug=True)
