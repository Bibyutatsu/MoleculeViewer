import csv
from flask import Flask, render_template, request, redirect, flash
from pager import Pager
from rdkit import Chem
from rdkit.Chem import Draw
import pandas as pd
import os
from moses.utils import get_mol
from pychem2 import constitution

APPNAME = "Molecule Viewer"
STATIC_FOLDER = 'example'
IMAGE_FOLDER = 'example/images'
TABLE_FILE = "example/fakecatalog.csv"

palette = ['#3fc5f0',
           '#42dee1',
           '#6decb9',
           '#eef5b2'
           ]


columns = ['name',
           'Weight',
           'AWeight',
           'nhyd',
           'nhal',
           'nhet',
           'nhev',
           'ncof',
           'ncocl',
           'ncobr',
           'ncoi',
           'ncarb',
           'nphos',
           'nsulph',
           'noxy',
           'nnitro',
           'nring',
           'nrot',
           'ndonr',
           'naccr',
           'nsb',
           'ndb',
           'naro',
           'ntb',
           'nta',
           'PC1',
           'PC2',
           'PC3',
           'PC4',
           'PC5',
           'PC6',
           'notes']


def read_table(url):
    """Return a list of dict"""
    # r = requests.get(url)
    with open(url) as f:
        return [row for row in csv.DictReader(f.readlines())]


def add_table(url, smiles):
    df = pd.read_csv(url)
    num_of_mols = len(df)
    mol = Chem.MolFromSmiles(smiles)
    notes = 'Molecule' + str(num_of_mols) + '.svg'
    out = constitution.GetConstitutional(mol)
    out['name'] = smiles
    out['notes'] = notes
    Draw.MolToFile(mol,
                   os.path.join(IMAGE_FOLDER, notes),
                   size=(1000, 1000), fitImage=True)
    row_to_add = pd.DataFrame([out])
    df = df.append(row_to_add)
    df.to_csv(url, index=False)


def delete_table():
    os.system('rm -r ./{0}/*'.format(IMAGE_FOLDER))
    pd.DataFrame([], columns=columns).to_csv(TABLE_FILE, index=False)


if os.path.isfile(TABLE_FILE):
    table = read_table(TABLE_FILE)
else:
    delete_table()
    table = read_table(TABLE_FILE)
pager = Pager(len(table))


app = Flask(__name__, static_folder=STATIC_FOLDER)
app.config.update(
    APPNAME=APPNAME,
)
app.secret_key = b'_5#y2L"F4Q8z\n\xec]/'


@app.after_request
def add_header(r):
    """
    Add headers to both force latest IE rendering engine or Chrome Frame,
    and also to cache the rendered page for 10 minutes.
    """
    r.headers["Cache-Control"] = "no-cache, no-store, must-revalidate"
    r.headers["Pragma"] = "no-cache"
    r.headers["Expires"] = "0"
    r.headers['Cache-Control'] = 'public, max-age=0'
    return r


@app.route('/')
def index():
    global pager, table
    return render_template('homeview.html', data=table)


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
            data=table[ind],
            column=columns,
            palette=palette)


@app.route('/goto', methods=['POST', 'GET'])
def goto():
    return redirect('/' + request.form['index'])


@app.route('/compareto', methods=['POST', 'GET'])
def compareto():
    return redirect(
        '/compare/' + request.form['index1'] + '&' + request.form['index2'])


@app.route('/molIn', methods=['POST', 'GET'])
def molIn():
    global table, pager
    if (
        (request.form['smiles'] != '') and
        (get_mol(request.form['smiles']) is not None)
    ):
        add_table(TABLE_FILE, request.form['smiles'])
        table = read_table(TABLE_FILE)
        pager = Pager(len(table))
    else:
        flash("Not a valid SMILE")
    return redirect('/')


@app.route('/molCsvIn', methods=['POST', 'GET'])
def molCsvIn():
    global table, pager
    if (
        (request.form['path'] != '') and
        (os.path.isfile(request.form['path']) is not False)
    ):
        df_inp = pd.read_csv(request.form['path'])
        for i, smile in enumerate(df_inp['SMILES']):
            if get_mol(smile) is not None:
                add_table(TABLE_FILE, smile)
            else:
                flash(str(i) + ' ' + smile + " is Not a valid SMILE")
            # percent = (i + 1) * 100 // len(df_inp)
            # render_template('progress.html', percent=percent)
        table = read_table(TABLE_FILE)
        pager = Pager(len(table))
    else:
        flash("Not a valid CSV path")
    return redirect('/')


@app.route('/clear')
def cleardb():
    global table, pager
    delete_table()
    table = read_table(TABLE_FILE)
    pager = Pager(len(table))
    return redirect('/')


@app.route('/compare/<int:ind1>&<int:ind2>/')
def compare(ind1, ind2):
    global table, pager
    if ind1 >= pager.count or ind2 >= pager.count:
        return render_template("404.html"), 404
    return render_template('compare.html',
                           pager=pager,
                           data1=table[ind1],
                           data2=table[ind2],
                           column=columns,
                           palette=palette)


if __name__ == '__main__':
    app.run(debug=True)
