import csv
from flask import Flask, render_template, request, redirect, flash
from pager import Pager
from rdkit import Chem
from rdkit.Chem import Draw
import pandas as pd
import os
from pychem2 import constitution, getmol
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
from radar import radar_factory

APPNAME = "Molecule Viewer"
STATIC_FOLDER = 'example'
MOLECULE_IMAGE_FOLDER = 'example/images'
RADAR_IMAGE_FOLDER = 'example/radars'
TABLE_FILE = "example/fakecatalog.csv"
COMPARISON = 'example/comparison'

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

radar_labels = ['Weight', 'AWeight', 'nhyd']


def create_radar(spoke_labels,
                 data,   # Should be a list of tuples with title and LOL data
                 labels,
                 save_path=os.path.join(RADAR_IMAGE_FOLDER,
                                        'Mol0.jpg'),
                 rgrids=[10, 50, 100, 500, 1000]):
    N = len(data[0][1][0])
    theta = radar_factory(N, frame='polygon')

    fig, axes = plt.subplots(figsize=(10, 10), nrows=len(data), ncols=1,
                             subplot_kw=dict(projection='radar'))
    fig.subplots_adjust(wspace=0.25, hspace=0.20, top=0.85, bottom=0.05)

    if len(data) == 1:
        axes = np.array([[axes]])

    colors = ['b', 'r', 'g', 'm', 'y']
    colors = colors[:len(data[0][1])]
    # Plot the four cases from the example data on separate axes
    for ax, (title, case_data) in zip(axes.flatten(), data):
        ax.set_rgrids(rgrids)
        ax.set_title(title, weight='bold', size='medium', position=(0.5, 1.1),
                     horizontalalignment='center', verticalalignment='center')
        for d, color in zip(case_data, colors):
            ax.plot(theta, d, color=color)
            ax.fill(theta, d, facecolor=color, alpha=0.25)
        ax.set_varlabels(spoke_labels)

    # add legend relative to top-left plot
    ax = axes[0, 0]
    ax.legend(labels, loc=(0.9, .95),
              labelspacing=0.1, fontsize='large')

    # fig.text(0.5, 0.965, '5-Factor Solution Profiles Across Four Scenarios',
    #          horizontalalignment='center', color='black', weight='bold',
    #          size='large')

    plt.savefig(save_path)


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

    data_radar = [[out[lab] for lab in radar_labels]]
    data = [(out['name'], data_radar)]

    create_radar(['Weight', 'Avg weight', 'NHYD'], data, ['base'],
                 save_path=os.path.join(RADAR_IMAGE_FOLDER,
                                        notes))
    Draw.MolToFile(mol,
                   os.path.join(MOLECULE_IMAGE_FOLDER, notes),
                   size=(1000, 1000), fitImage=True, imageType='svg')
    row_to_add = pd.DataFrame([out])
    df = df.append(row_to_add)
    df.to_csv(url, index=False)


def delete_table():
    os.system('rm -r ./{0}/*'.format(MOLECULE_IMAGE_FOLDER))
    os.system('rm -r ./{0}/*'.format(RADAR_IMAGE_FOLDER))
    os.system('rm -r ./{0}/*'.format(COMPARISON))
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
        (getmol.ReadMolFromSmile(request.form['smiles']) is not None)
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
            if getmol.ReadMolFromSmile(smile) is not None:
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

    data_radar = [[table[inds][lab] for lab in radar_labels
                   ] for inds in [ind1, ind2]]
    data = [('base', data_radar)]
    create_radar(['Weight', 'Avg weight', 'NHYD'], data,
                 [table[ind1]['name'], table[ind2]['name']],
                 save_path=os.path.join(COMPARISON,
                                        str(ind1) + '&' + str(ind2) + '.svg'))

    if ind1 >= pager.count or ind2 >= pager.count:
        return render_template("404.html"), 404
    return render_template('compare.html',
                           pager=pager,
                           data1=table[ind1],
                           data2=table[ind2],
                           column=columns,
                           palette=palette,
                           ind1=str(ind1),
                           ind2=str(ind2))


if __name__ == '__main__':
    app.run(debug=True)
