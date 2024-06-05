// script.js

document.addEventListener('DOMContentLoaded', function() {
    document.querySelectorAll('.uploadForm').forEach(form => {
        form.addEventListener('submit', onFormSubmit);
    });
    document.getElementById('uploadScatteringForm').addEventListener('submit', onScatteringFormSubmit);
    console.log("Event listeners added to forms");
});

function onFormSubmit(event) {
    event.preventDefault();
    console.log("Form submitted");
    handleFileUpload(event.target);
}

function onScatteringFormSubmit(event) {
    event.preventDefault();
    console.log("Scattering form submitted");
    handleScatteringFileUpload(event.target);
}

function handleFileUpload(form) {
    var fileInput = form.querySelector('.fileInput');
    var file = fileInput.files[0];
    if (!file) {
        alert("Please select a file.");
        return;
    }

    var chainType = form.querySelector('input[name="chainType"]:checked').value;

    var formData = new FormData();
    formData.append('file', file);
    formData.append('chainType', chainType);

    var xhr = new XMLHttpRequest();
    xhr.open('POST', '/upload', true);

    xhr.onload = function () {
        if (xhr.status === 200) {
            var response = JSON.parse(xhr.responseText);
            if (response.chainType === 'monomer') {
                response.filenames.forEach(function(filename) {
                    console.log("Chain file uploaded successfully:", filename);
                    createNGLViewer(filename, file.name, filename.split('_').pop().split('.')[0]);
                });
            } else if (response.chainType === 'multimer') {
                console.log("File uploaded successfully:", response.filename);
                createNGLViewer(response.filename, file.name);
            }
            addNewUploadForm();
            form.querySelector('input[type="submit"]').disabled = true; // Disable the submit button to prevent resubmission
        } else {
            alert('An error occurred while uploading the file.');
        }
    };

    xhr.onerror = function () {
        alert('An error occurred while uploading the file.');
    };

    xhr.send(formData);
}

function handleScatteringFileUpload(form) {
    var fileInput = form.querySelector('.scatteringFileInput');
    var file = fileInput.files[0];
    if (!file) {
        alert("Please select a scattering profile file.");
        return;
    }

    var formData = new FormData();
    formData.append('scatteringFile', file);

    var xhr = new XMLHttpRequest();
    xhr.open('POST', '/uploadScattering', true);

    xhr.onload = function () {
        if (xhr.status === 200) {
            var response = JSON.parse(xhr.responseText);
            if (response.data) {
                console.log("Scattering file uploaded successfully");
                plotScatteringProfile(response.data);
            } else {
                alert("Error: " + response.error);
            }
        } else {
            alert('An error occurred while uploading the scattering profile.');
        }
    };

    xhr.onerror = function () {
        alert('An error occurred while uploading the scattering profile.');
    };

    xhr.send(formData);
}

let currentStage = null;

function createNGLViewer(filename, originalFilename, chainId) {
    var container = document.createElement('div');
    container.className = 'viewport';
    container.id = 'viewport-' + filename;

    var title = document.createElement('div');
    title.className = 'viewer-title';
    title.textContent = chainId ? `${originalFilename} - Chain ${chainId}` : originalFilename;

    var wrapper = document.createElement('div');
    wrapper.appendChild(title);
    wrapper.appendChild(container);

    document.getElementById('viewports').appendChild(wrapper);

    currentStage = new NGL.Stage(container, { backgroundColor: 'white' });
    currentStage.loadFile("/uploads/" + filename).then(function(o) {
        console.log("File loaded into NGL:", filename);
        applyColorScheme(o);
        o.autoView();
    }).catch(function(error) {
        console.error("Error loading file into NGL:", error);
        alert('An error occurred while loading the file into NGL.');
    });
}

function applyColorScheme(object) {
    let colorDropdown = document.getElementById('colorDropdown');
    let selectedOption = colorDropdown.value;

    object.addRepresentation("cartoon", {
        color: function(atom) {
            if (selectedOption === 'helix' && atom.sstruc === 'h') {
                return '#ff0000'; // red
            } else if (selectedOption === 'sheet' && atom.sstruc === 's') {
                return '#ff0000'; // red
            } else if (selectedOption === 'coil' && !atom.sstruc) {
                return '#ff0000'; // red
            } else {
                return '#aaaaaa'; // grey
            }
        }
    });
}

function updateColorScheme() {
    if (currentStage) {
        currentStage.removeAllComponents();
        currentStage.eachComponent(function(component) {
            applyColorScheme(component);
        });
    }
}

function addNewUploadForm() {
    var newFormId = 'uploadForm-' + (document.querySelectorAll('.uploadForm').length + 1);
    var newForm = document.createElement('form');
    newForm.id = newFormId;
    newForm.className = 'uploadForm';
    newForm.enctype = 'multipart/form-data';
    newForm.innerHTML = `
        <input type="file" name="file" class="fileInput">
        <div>
            <input type="radio" name="chainType" value="monomer" checked> Monomer
            <input type="radio" name="chainType" value="multimer"> Multimer
        </div>
        <input type="submit" value="Upload">
    `;
    newForm.addEventListener('submit', onFormSubmit);
    document.getElementById('uploadForms').appendChild(newForm);
    console.log("New upload form added");
}

function plotScatteringProfile(data) {
    var trace = {
        x: data.q,
        y: data.intensity,
        mode: 'markers',  // scatter plot with markers
        type: 'scatter'
    };

    var layout = {
        title: 'Scattering Profile',
        xaxis: {
            title: 'q'
        },
        yaxis: {
            title: 'Intensity'
        },
        width: 600,  // limit the width of the plot
        height: 400  // limit the height of the plot
    };

    Plotly.newPlot('scatteringPlot', [trace], layout);
}
