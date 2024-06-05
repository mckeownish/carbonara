from flask import Flask, request, session, render_template, jsonify, send_from_directory
from werkzeug.utils import secure_filename
import os
from pdbfixer import PDBFixer
from openmm.app import PDBFile
import numpy as np

app = Flask(__name__)
app.config['UPLOAD_FOLDER'] = 'uploads'
app.config['SECRET_KEY'] = 'supersecretkey'
app.config['SESSION_TYPE'] = 'filesystem'

if not os.path.exists(app.config['UPLOAD_FOLDER']):
    os.makedirs(app.config['UPLOAD_FOLDER'])

def allowed_file(filename):
    return '.' in filename and filename.rsplit('.', 1)[1].lower() in {'pdb', 'dat'}

def is_valid_pdb(pdb_path):
    try:
        pdb = PDBFile(pdb_path)
        return True
    except Exception as e:
        print(f"Error reading PDB file: {e}")
        return False

def split_pdb_into_chains(pdb_file_path):
    from Bio import PDB
    parser = PDB.PDBParser()
    structure = parser.get_structure('structure', pdb_file_path)
    chain_files = []

    for model in structure:
        for chain in model:
            io = PDB.PDBIO()
            io.set_structure(chain)
            chain_file_path = pdb_file_path.replace('.pdb', f'_{chain.id}.pdb')
            io.save(chain_file_path)
            chain_files.append(chain_file_path)

    return chain_files

@app.route('/')
def index():
    session.clear()
    return render_template('upload.html')

@app.route('/upload_pdb', methods=['POST'])
def upload_pdb():
    if 'pdbFile' not in request.files:
        return jsonify({'error': 'No PDB file part'}), 400
    
    pdb_file = request.files['pdbFile']
    print(request.files)
    
    if pdb_file.filename == '':
        return jsonify({'error': 'No selected PDB file'}), 400

    pdb_filename = secure_filename(pdb_file.filename)
    pdb_path = os.path.join(app.config['UPLOAD_FOLDER'], pdb_filename)
    print("******" + pdb_path)
    pdb_file.save(pdb_path)
    
    if not is_valid_pdb(pdb_path):
        return jsonify({'error': 'Invalid PDB file'}), 400

    chain_type = request.form.get('chainType')
    
    new_files = []
    if chain_type == 'monomer':
        try:
            chain_files = split_pdb_into_chains(pdb_path)
            new_files.extend(chain_files)
        except Exception as e:
            return jsonify({'error': f'Error splitting PDB file into chains: {e}'}), 500
    else:
        new_files.append(pdb_path)
    
    if 'pdb_files' not in session:
        session['pdb_files'] = []
    if 'no_structures' not in session:
        session['no_structures'] = 0

    session['pdb_files'].extend(new_files)
    session['no_structures'] += len(new_files)
    
    return jsonify({'pdb_files': session['pdb_files']}), 200

@app.route('/uploads/<filename>')
def uploaded_file(filename):
    return send_from_directory(app.config['UPLOAD_FOLDER'], filename)

@app.route('/upload_saxs', methods=['POST'])
def upload_saxs():
    if 'scatteringFile' not in request.files:
        return jsonify({'error': 'No file part'}), 400

    file = request.files['scatteringFile']

    if file.filename == '':
        return jsonify({'error': 'No selected file'}), 400

    if file and allowed_file(file.filename):
        filename = file.filename
        file_path = os.path.join(app.config['UPLOAD_FOLDER'], filename)
        file.save(file_path)

        processed_file_path = write_saxs(file_path, app.config['UPLOAD_FOLDER'])
        data = pd.read_csv(processed_file_path, delim_whitespace=True, names=['q', 'intensity'])
        data_dict = data.to_dict(orient='list')
        return jsonify({'data': data_dict}), 200

    return jsonify({'error': 'File not allowed'}), 400



def write_saxs(SAXS_file, working_path):
    with open(SAXS_file) as oldfile, open('temp.txt', 'w') as newfile:
        for line in oldfile:
            final_list = []
            for elem in line.split():
                try:
                    float(elem)
                except ValueError:
                    final_list.append(elem)
            if len(final_list) == 0:
                newfile.write(line)

    saxs_arr = np.genfromtxt('temp.txt')

    if saxs_arr.shape[1] == 3:
        saxs_arr = saxs_arr[:, :2]

    # Check if it is in angstroms; if the last value is >1 we assume it's in nanometers.
    if saxs_arr[-1, 0] > 1:
        for i in range(0, len(saxs_arr)):
            saxs_arr[i, 0] = saxs_arr[i, 0] / 10.0

    file_path_name = os.path.join(working_path, 'Saxs.dat')
    np.savetxt(file_path_name, saxs_arr, delimiter=' ', fmt='%s', newline='\n', header='', footer='')
    os.remove("temp.txt")

    return file_path_name



@app.route('/configure', methods=['GET', 'POST'])
def configure_parameters():
    if request.method == 'POST':
        session['save_location'] = request.form['save_location']
        session['constrained_distances'] = request.form['constrained_distances']
        session['qmin'] = request.form['qmin']
        session['qmax'] = request.form['qmax']
        session['fit_steps'] = request.form['fit_steps']
        return redirect(url_for('generate_script'))
    
    return render_template('configure.html', no_structures=session['no_structures'])

@app.route('/generate_script')
def generate_script():
    save_location = session.get('save_location')
    pdb_files = session.get('pdb_files')
    saxs_file = session.get('saxs_file')
    constrained_distances = session.get('constrained_distances')
    qmin = session.get('qmin')
    qmax = session.get('qmax')
    fit_steps = session.get('fit_steps')
    no_structures = session.get('no_structures')

    script_content = f"""
#!/bin/bash

ROOT=$(dirname "$(readlink -f "$0")")

ScatterFile={saxs_file}
fileLocs={','.join(pdb_files)}
initialCoordsFile=frompdb
pairedPredictions=False
fixedsections={constrained_distances}
noStructures={no_structures}
withinMonomerHydroCover=none
betweenMonomerHydroCover=none
kmin={qmin}
kmax={qmax}
maxNoFitSteps={fit_steps}
predictionFile={save_location}/fitdata
scatterOut={save_location}/fitdata
mixtureFile={save_location}/mixtureFile.dat
prevFitStr={save_location}/redundant
logLoc={save_location}/fitdata
endLinePrevLog=null
affineTrans=False

for i in {{1..{no_structures}}}
do
    echo "\\n"
    echo " >> Run number : $i "
    echo "\\n"
    ./predictStructureQvary $ScatterFile $fileLocs $initialCoordsFile $pairedPredictions $fixedsections $noStructures $withinMonomerHydroCover $betweenMonomerHydroCover $kmin $kmax $maxNoFitSteps $predictionFile/mol$i $scatterOut/scatter$i.dat $mixtureFile $prevFitStr $logLoc/fitLog$i.dat $endLinePrevLog $affineTrans
done
"""

    script_path = os.path.join(app.config['UPLOAD_FOLDER'], 'run_script.sh')
    with open(script_path, 'w') as script_file:
        script_file.write(script_content)

    return send_file(script_path, as_attachment=True, download_name='run_script.sh')

if __name__ == '__main__':
    app.run(debug=True, port=5001)  # Changed port to avoid conflict
