/* Implicit hydration WAXSiS approximation */

#include "ktlMoleculeRandom.h"
#include <vector>
#include <map>
#include <string>
#include <sstream>

#include <cmath>
#include "writheFP.h"


struct ScatteringCenters {
    std::vector<point> coordinates;
    std::vector<char> types;
    std::vector<std::vector<double>> distances;
};

struct ExperimentalData {
   std::vector<double> q;
   std::vector<double> I;
   std::vector<double> I_err;
};



void readInStructures(std::vector<ktlMolecule>& mol, int noStructures, std::string path) {

    if (noStructures <= 0) {
        std::cerr << "Invalid number of structures specified." << std::endl;
        return;
    }
    
    for (int i = 0; i < noStructures; i++) {
        
        ktlMolecule molTmp;
        
        double rmin = 3.7;
        double rmax = 3.9;
        double lmin = 4.0;
        
        // Is this right? fingerPrint1 should contain all chain info
        std::string sequenceLoc = path + "fingerPrint" + std::to_string(i + 1) + ".dat";
        molTmp.readInSequence(sequenceLoc.c_str(), rmin, rmax,  lmin);
        
        std::string coordinateLoc = path + "coordinates" + std::to_string(i + 1) + ".dat";
        molTmp.readInCoordinates(coordinateLoc.c_str());
        
        molTmp.getHydrophobicResidues();
        
        mol.push_back(molTmp);
    }
}


std::vector<char> extract_sequence(const char* filename) {
   std::ifstream file(filename);
   std::string line;
   
   if (!file.is_open()) {
       std::cerr << "Could not open file: " << filename << std::endl;
       return std::vector<char>();
   }
   
   std::getline(file, line); // Skip first line
   std::getline(file, line); // Skip empty line
   std::getline(file, line); // Get sequence
   
   return std::vector<char>(line.begin(), line.end());
}


std::vector<point> flatten_coords( std::vector<std::vector<point>>& coords) {
   std::vector<point> flat;
   for (const auto& section : coords) {
       flat.insert(flat.end(), section.begin(), section.end());
   }
   return flat;
}


std::vector<point> calculate_geometric_normals(std::vector<point>& ca_coords) {
    
    std::vector<point> ca_vectors;
    std::vector<point> normals;
    
    // Calculate vectors between consecutive CA atoms
    for (size_t i = 1; i < ca_coords.size(); i++) {
        ca_vectors.push_back(ca_coords[i] - ca_coords[i-1]);
    }
    
    // calculate the geometric *outward* facing normals
    for (size_t i = 1; i < ca_vectors.size(); i++) {
        
        point norm = ca_vectors[i] - ca_vectors[i-1];
        norm = norm * (-1.0) / norm.length();   // outward
        normals.push_back(norm);
        
    }
    
    point first_norm = ca_vectors[0] * (-1.0) / ca_vectors[0].length();
    point final_norm = ca_vectors.back() / ca_vectors.back().length();
    
    normals.insert( normals.begin(), first_norm );
    normals.push_back( final_norm );
    
    return normals;
}


std::vector<point> place_side_chains( std::vector<point>& ca_coords,
                                      std::vector<point>& geometric_vectors,
                                      std::vector<char>& residue_names) {
    
    std::map<char, double> residue_distances = {
        {'R', 4.2662}, {'N', 2.5349}, {'D', 2.5558}, {'C', 2.3839},
        {'Q', 3.1861}, {'E', 3.2541}, {'H', 3.1861}, {'I', 2.3115},
        {'L', 2.6183}, {'K', 3.6349}, {'M', 3.1912}, {'F', 3.4033},
        {'P', 1.8773}, {'U', 1.5419}, {'S', 1.9661}, {'T', 1.9533},
        {'W', 3.8916}, {'Y', 3.8807}, {'V', 1.9555}, {'G', 0.0}, {'A', 0.0}
    };
    
    std::vector<point> side_chain_positions;
    
    for (size_t i = 0; i < ca_coords.size(); i++) {
        
        double distance = residue_distances[residue_names[i]];
        side_chain_positions.push_back(ca_coords[i] + geometric_vectors[i] * distance);
        
    }
    
    return side_chain_positions;
}


std::vector<std::vector<double>> calculate_distances( std::vector<point>& coordinates ) {
    
    size_t n = coordinates.size();
    std::vector< std::vector<double> > distances(n, std::vector<double>(n));
    double min_dist = std::numeric_limits<double>::max();
    
    for (size_t i = 0; i < n; i++){
        for (size_t j = i; j < n; j++){
            
            if (i==j) { distances[i][j] = 0.0;
                
            } else {
                double dist = coordinates[i].eDist(coordinates[j]);
                distances[i][j] = dist;
                distances[j][i] = dist;
            }
        }
    }
    
    return distances;
}

double find_minimum_distance( std::vector< std::vector<double> > distances ) {

    double min_dist = std::numeric_limits<double>::max();

    for (size_t i = 0; i < distances.size(); i++){
        for (size_t j = i; j < distances.size(); j++){

            if (i!=j) { if (distances[i][j] < min_dist) { min_dist = distances[i][j]; } }

        }
    }

    return min_dist;
}


ScatteringCenters process_structure(std::vector<point>& ca_coords,
                                    std::vector<char>& residue_names){
    
    std::vector<point> geometric_normals = calculate_geometric_normals(ca_coords);
    std::vector<point> side_chain_positions = place_side_chains(ca_coords, geometric_normals, residue_names);
        
    ScatteringCenters centers;
    
    // add backbone + sidechains centres for non-GLY/ALA residues
    for (size_t i=0; i < residue_names.size(); i++) {
        
        if (residue_names[i] != 'A' && residue_names[i] != 'G') {
            
            // backbone
            centers.coordinates.push_back(ca_coords[i]);
            centers.types.push_back('B');
            
            // side chain
            centers.coordinates.push_back(side_chain_positions[i]);
            centers.types.push_back(residue_names[i]);
        }
    }
    
    // add GLY + ALA CA postions
    for (size_t i=0; i < residue_names.size(); i++) {
        
        if (residue_names[i] == 'A' || residue_names[i] == 'G') {
            
            // backbone
            centers.coordinates.push_back(ca_coords[i]);
            centers.types.push_back(residue_names[i]);
        }
    }
    
    centers.distances = calculate_distances(centers.coordinates);
    
    return centers;
}


std::map< char, std::vector<double> > init_form_factors() {
    
    std::map< char, std::vector<double> > form_factors;
    
    form_factors['B'] = {5.325233, 5.321616, 5.321073, 5.330379, 5.367118, 5.418221, 5.548391, 5.694481, 5.802268, 5.880531, 5.900372, 5.865756, 5.778460, 5.611933, 5.353311, 5.065751, 4.844303, 4.655859, 4.507514, 4.382296, 4.200807};
    form_factors['A'] = {6.860783, 6.884742, 6.936093, 6.987066, 7.002729, 7.001643, 7.011408, 7.059012, 7.136653, 7.204735, 7.278987, 7.358991, 7.468390, 7.587982, 7.698403, 7.792850, 7.835658, 7.747937, 7.795444, 8.008073, 8.268595};
    form_factors['R'] = {15.586594, 15.603896, 15.644747, 15.706051, 15.731217, 15.786658, 15.751397, 15.723773, 15.667085, 15.568371, 15.501597, 15.510941, 15.660628, 15.938443, 16.391161, 17.138668, 18.150400, 19.201813, 19.957359, 20.556793, 21.194653};
    form_factors['N'] = {14.042562, 14.039913, 14.030098, 14.046182, 14.072546, 14.079157, 13.947478, 13.831409, 13.740083, 13.667616, 13.673161, 13.796733, 14.041808, 14.428295, 14.873251, 15.052513, 14.851926, 14.920277, 15.730271, 17.125431, 18.337770};
    form_factors['D'] = {18.276451, 18.266953, 18.238106, 18.188072, 18.080282, 17.894873, 17.509510, 17.006626, 16.529716, 16.041245, 15.567556, 15.144760, 14.858373, 14.703475, 14.630036, 14.854190, 15.309014, 15.771427, 16.001890, 16.368433, 16.704248};
    form_factors['C'] = {10.351711, 10.346272, 10.323308, 10.293743, 10.243892, 10.212900, 10.157944, 10.163668, 10.267132, 10.432532, 10.664062, 10.914792, 11.118233, 11.204067, 11.197291, 11.052269, 10.771274, 10.325391, 9.898160, 9.588455, 9.458176};
    form_factors['Q'] = {13.376416, 13.384198, 13.387307, 13.362279, 13.276151, 13.152504, 12.912460, 12.665339, 12.532801, 12.473034, 12.546851, 12.777804, 13.243794, 13.995014, 15.000033, 15.920972, 16.500341, 17.029819, 17.050764, 16.716652, 16.530117};
    form_factors['E'] = {21.686304, 21.668835, 21.592291, 21.449541, 21.200043, 20.835613, 20.150393, 19.448851, 18.805557, 18.179039, 17.611521, 17.140610, 16.856714, 16.806103, 16.947565, 17.285227, 17.668753, 18.353695, 19.242361, 19.913992, 20.734098};
    form_factors['G'] = {10.030974, 10.004344, 9.943201, 9.875965, 9.802672, 9.715462, 9.611974, 9.495539, 9.422803, 9.374266, 9.338644, 9.343287, 9.433417, 9.612279, 9.818403, 10.084248, 10.328732, 10.578753, 10.550336, 10.467358, 10.294571};
    form_factors['H'] = {13.493807, 13.499746, 13.497142, 13.453699, 13.357458, 13.200852, 12.953443, 12.703934, 12.458268, 12.168073, 11.901849, 11.685180, 11.584588, 11.724971, 12.122943, 12.761358, 13.564467, 14.011192, 14.246178, 14.511696, 15.254260};
    form_factors['I'] = {0.027768, 0.028121, 0.028573, 0.023142, 0.007413, 0.003274, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.010318, 0.029916, 0.021361, 0.011234, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000};
    form_factors['L'] = {0.001379, 0.001181, 0.000961, 0.000118, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.012047, 0.007686, 0.000000, 0.000000, 0.001314, 0.006526, 0.009532, 0.008636};
    form_factors['K'] = {9.349114, 9.359500, 9.370693, 9.372478, 9.334329, 9.252807, 9.074130, 8.906091, 8.755677, 8.595029, 8.477226, 8.410988, 8.422702, 8.465599, 8.387974, 8.202576, 7.873840, 7.777987, 8.014545, 8.801325, 10.069035};
    form_factors['M'] = {5.460539, 5.465856, 5.488052, 5.560043, 5.633809, 5.795099, 6.002742, 6.323314, 6.623469, 6.956409, 7.325475, 7.701069, 8.105552, 8.563451, 9.163991, 9.708839, 9.891646, 9.833041, 9.734925, 9.476597, 9.286131};
    form_factors['F'] = {0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.001274, 0.020930, 0.086051, 0.251851, 0.409124, 0.562665, 0.686578, 0.832251, 1.048074, 1.359721, 1.720046, 2.006726, 2.087998, 2.053319, 2.235860, 2.837866};
    form_factors['P'] = {4.716053, 4.711367, 4.685131, 4.642177, 4.559615, 4.457505, 4.098892, 3.705757, 3.372101, 3.053588, 2.797038, 2.656890, 2.760468, 3.205317, 3.912106, 4.736746, 5.429208, 6.331475, 7.254241, 7.861680, 7.935802};
    form_factors['S'] = {8.985313, 9.003293, 9.028045, 9.040872, 8.994022, 8.931241, 8.708052, 8.462237, 8.256139, 8.054568, 7.901446, 7.829412, 7.892257, 8.112808, 8.520312, 9.034179, 9.444173, 9.954925, 10.321150, 9.990824, 9.688856};
    form_factors['T'] = {6.628838, 6.637546, 6.641963, 6.629081, 6.596369, 6.589334, 6.517858, 6.426005, 6.370910, 6.309315, 6.288209, 6.310956, 6.417308, 6.626843, 6.961301, 7.252404, 7.440239, 7.624983, 7.722063, 7.757763, 8.270276};
    form_factors['W'] = {2.888272, 2.949734, 3.104849, 3.325679, 3.557907, 3.812995, 4.198963, 4.642990, 5.063122, 5.487359, 5.964144, 6.477489, 7.032223, 7.543501, 7.751148, 7.804211, 7.761533, 8.023214, 8.654520, 9.687876, 10.764302};
    form_factors['Y'] = {5.100512, 5.130500, 5.228212, 5.390429, 5.591930, 5.801286, 5.978783, 6.084779, 6.201869, 6.252052, 6.266182, 6.257378, 6.282480, 6.352914, 6.603123, 6.989089, 7.407878, 7.660326, 7.929643, 8.116093, 8.308664};
    form_factors['V'] = {0.028584, 0.028193, 0.026019, 0.023774, 0.010495, 0.002455, 0.000000, 0.000000, 0.006006, 0.012980, 0.016876, 0.014590, 0.028809, 0.029842, 0.008057, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000};
    
    return form_factors;
}


std::vector<double> init_q_values() {
    
    std::vector<double> q;
    for (int i = 0; i <= 20; i++) {
        q.push_back(0.01 * i);
    }
    
    return q;
}


std::vector<double> calculate_saxs_implicit( ScatteringCenters centers ) {
    
    auto q_values = init_q_values();
    auto form_factors = init_form_factors();
    
    std::vector<double> I(q_values.size(), 0.0);
    
    std::vector< std::vector<double> > ff_matrix;
    for (char type : centers.types) {
        ff_matrix.push_back( form_factors.at(type) );
    }
    
    // calculate scattering at each q
    for (size_t q_idx = 0; q_idx < q_values.size(); q_idx++) {
        
        double q = q_values[q_idx];
        
        for (size_t i = 0; i < centers.coordinates.size(); i++) {
            for (size_t j = 0; j < centers.coordinates.size(); j++) {
                double r_ij = centers.distances[i][j];
                double ff_product = ff_matrix[i][q_idx] * ff_matrix[j][q_idx];
                
                double sinc;
                if (q == 0) {
                    sinc = 1.0;
                } else {
                    sinc = (r_ij == 0) ? 1.0 : std::sin(q * r_ij) / (q * r_ij);
                }
                
                I[q_idx] += ff_product * sinc;
            }
        }
    }
    
    return I;
}


void calculate_spline_coefficients(const std::vector<double>& x, const std::vector<double>& y,
                                 std::vector<double>& A, std::vector<double>& B,
                                 std::vector<double>& C, std::vector<double>& D) {
    size_t n = x.size();
    std::vector<double> h(n - 1), alpha(n - 1), l(n), mu(n), z(n);
    
    // Step 1: Compute h_i = x_{i+1} - x_i
    for (size_t i = 0; i < n - 1; i++) {
        h[i] = x[i + 1] - x[i];
    }
    
    // Step 2: Compute alpha_i
    for (size_t i = 1; i < n - 1; i++) {
        alpha[i] = (3.0/h[i])*(y[i+1] - y[i]) - (3.0/h[i-1])*(y[i] - y[i-1]);
    }
    
    // Step 3: Solve tridiagonal system
    l[0] = 1.0;
    mu[0] = 0.0;
    z[0] = 0.0;
    
    for (size_t i = 1; i < n - 1; i++) {
        l[i] = 2.0 * (x[i+1] - x[i-1]) - h[i-1]*mu[i-1];
        mu[i] = h[i]/l[i];
        z[i] = (alpha[i] - h[i-1]*z[i-1])/l[i];
    }
    
    l[n-1] = 1.0;
    z[n-1] = 0.0;
    
    // Step 4: Compute coefficients
    std::vector<double> c(n);
    c[n-1] = 0.0;
    
    for (size_t j = n - 2; j < n; j--) {  // Note: using size_t, so we break when j underflows
        c[j] = z[j] - mu[j]*c[j+1];
    }
    
    // Step 5: Compute spline coefficients
    for (size_t i = 0; i < n - 1; i++) {
        A[i] = (c[i+1] - c[i])/(3.0*h[i]);
        B[i] = c[i];
        C[i] = (y[i+1] - y[i])/h[i] - h[i]*(c[i+1] + 2.0*c[i])/3.0;
        D[i] = y[i];
    }
}

double evaluate_spline(double x_val, const std::vector<double>& x,
                      const std::vector<double>& A, const std::vector<double>& B,
                      const std::vector<double>& C, const std::vector<double>& D) {
    // Find the appropriate interval
    size_t i = 0;
    while (i < x.size() - 2 && x[i + 1] <= x_val) i++;
    
    // Compute relative x
    double dx = x_val - x[i];
    
    // Evaluate cubic polynomial
    return ((A[i]*dx + B[i])*dx + C[i])*dx + D[i];
}


std::vector<double> calculate_intensity_at_experimental_q(  const ExperimentalData& exp_data,
                                                            const std::vector<double>& q_mod,
                                                            const std::vector<double>& I_mod) {
    
    double q_max = std::min(0.2, *std::max_element(exp_data.q.begin(), exp_data.q.end()));
    
    size_t N = q_mod.size() - 1;
    std::vector<double> A(N), B(N), C(N), D(N);
    
    calculate_spline_coefficients(q_mod, I_mod, A, B, C, D);
    
    std::vector<double> I_interp;
    for(double q : exp_data.q) {
        
        if (q <= q_max) {
            I_interp.push_back(evaluate_spline(q, q_mod, A, B, C, D));
        }
        
    }
    return I_interp;
}


ExperimentalData read_experimental_data(const char* filename) {
   ExperimentalData data;
   std::ifstream file(filename);
   double q, I, I_err;
   
   while (file >> q >> I >> I_err) {
       
       if (q <= 0.2) {
           data.q.push_back(q);
           data.I.push_back(I);
           data.I_err.push_back(I_err);
       }
   }
   return data;
}


double calculate_scale_factor( const ExperimentalData& exp_data,
                             const std::vector<double>& I_mod_interp,
                             bool use_errors = false ) {
                                
    //scale = Σ( Iexp * Imod / σ² ) / Σ( Imod² / σ² )  [if use_errors]
    //scale = Σ( Iexp * Imod ) / Σ( Imod² )            [if !use_errors]

    double numerator = 0.0;
    double denominator = 0.0;

    if (exp_data.q.size() != I_mod_interp.size()) {
        std::cerr << "Error: Experimental data and interpolated model data sizes do not match. ";
        std::cerr << "exp_data.q size: " << exp_data.q.size() << ", ";
        std::cerr << "I_mod_interp size: " << I_mod_interp.size() << std::endl;
        return 0.0;
    }
    
    for (size_t i = 0; i < exp_data.q.size(); ++i) {
        double weight = use_errors ? 1.0 / (exp_data.I_err[i] * exp_data.I_err[i]) : 1.0;
        numerator += exp_data.I[i] * I_mod_interp[i] * weight;
        denominator += I_mod_interp[i] * I_mod_interp[i] * weight;
    }

    return numerator / denominator;
}


double calculate_chi_squared( const ExperimentalData exp_data,
                              const std::vector<double>& I_mod_interp,
                              double scale_factor ) {
    
    double chi_squared = 0.0;

    for (size_t i = 0; i < exp_data.q.size(); ++i) {

        double diff = exp_data.I[i] - scale_factor * I_mod_interp[i];
        chi_squared += diff * diff / (exp_data.I_err[i] * exp_data.I_err[i]);
    }

    return chi_squared / (exp_data.q.size() - 1);
}


void write_to_file( const std::string& filename,
                    const ExperimentalData exp_data,
                    const std::vector<double>& I_mod_interp,
                    double scale_factor,
                    double chi_squared ) {

    std::ofstream outfile(filename);

    if (!outfile.is_open(   )) {
        std::cerr << "Failed to open output file" << std::endl;
        return;
    }

    outfile << "Chi squared: " << chi_squared << std::endl;
    outfile << "Scale factor: " << scale_factor << std::endl;

    outfile << "q, I_exp, I_err, I_mod_scaled" << std::endl;

    for (size_t i = 0; i < exp_data.q.size(); ++i) {
        outfile << exp_data.q[i] << " "
                << exp_data.I[i] << " "
                << exp_data.I_err[i] << " "
                << scale_factor * I_mod_interp[i] << std::endl;
    }

    outfile.close();
    
}


std::vector<int> readFlexibleSections(const std::string& path, int structureIndex) {
    std::vector<int> flexible_sections;
    std::string filename = path + "varyingSectionSecondary" + std::to_string(structureIndex + 1) + ".dat";
    
    std::ifstream file(filename);
    if (!file.is_open()) {
        std::cerr << "Failed to open varying section file: " << filename << std::endl;
        return flexible_sections;
    }

    std::string line;
    while (std::getline(file, line) && !line.empty()) {
        std::stringstream ss(line);
        int index;
        while (ss >> index) {
            flexible_sections.push_back(index);
        }
    }
    
    file.close();
    return flexible_sections;
}


struct SectionInfo {
    int chain_number;
    int section_index;
    int global_index;
};

std::vector<SectionInfo> get_flexible_indices( ktlMolecule& mol, 
                                               std::vector<int>& vary_sections,
                                               bool vary_all = false ) {
    std::vector<SectionInfo> sections;
    int global_index = 0;
    
    for (int chain = 1; chain <= mol.noChains(); chain++) {
        for (int section = 0; section < mol.getSubsecSize(chain) - 1; section++) {
            if (vary_all || std::find(vary_sections.begin(), vary_sections.end(), global_index) != vary_sections.end()) {
                sections.push_back({chain, section, global_index});
            }
            global_index++;
        }
    }

    return sections;
}

void updateProgressBar(int current, int total, int barWidth = 50) {
    float progress = static_cast<float>(current) / total;
    int pos = static_cast<int>(barWidth * progress);
    
    std::cout << "\r[";
    for (int i = 0; i < barWidth; ++i) {
        if (i < pos) std::cout << "=";
        else if (i == pos) std::cout << ">";
        else std::cout << " ";
    }
    std::cout << "] " << int(progress * 100.0) << "% (" << current << "/" << total << ")" << std::flush;
}


/* Writhe stuff */


// Helper function to calculate unit vector
point unit_vec( point& v1, point& v2) {
    point cross = v1.cross(v2);
    double length = cross.length();
    if (length < 1e-10) {
        return point(0, 0, 0);
    }
    return cross * (1.0/length);
}

// Calculate writhe contribution between 4 points
double Gauss_Int_4_segment( point& p1,  point& p2, point& p3, point& p4) {
    // Calculate vectors between points
    point r12 = p2 - p1;
    point r13 = p3 - p1;
    point r14 = p4 - p1;
    point r23 = p3 - p2;
    point r24 = p4 - p2;
    point r34 = p4 - p3;

    // Calculate normal vectors (unit vectors perpendicular to triangle faces)
    point n1 = unit_vec(r13, r14);  // Normal to triangle p1-p3-p4
    point n2 = unit_vec(r14, r24);  // Normal to triangle p1-p4-p2
    point n3 = unit_vec(r24, r23);  // Normal to triangle p2-p4-p3
    point n4 = unit_vec(r23, r13);  // Normal to triangle p2-p3-p1

    // Skip if any normal vectors are zero
    if (!n1.isNonzero() || !n2.isNonzero() || !n3.isNonzero() || !n4.isNonzero()) {
        return 0.0;
    }

    // Calculate solid angle using asin of dot products
    double Sigma_star = asin(n1.dotprod(n2)) + asin(n2.dotprod(n3)) 
                     + asin(n3.dotprod(n4)) + asin(n4.dotprod(n1));
    
    // Determine sign based on orientation
    double sign = (r34.cross(r12)).dotprod(r13) > 0 ? 1.0 : -1.0;
    
    return (1.0/(4.0*M_PI)) * Sigma_star * sign;
}

// Function to calculate the full writhe matrix
std::vector<std::vector<double>> calculateWritheMatrix( std::vector<point>& coords) {
    int n_segments = coords.size() - 1;
    std::vector<std::vector<double>> writhe_matrix(n_segments, 
                                                  std::vector<double>(n_segments, 0.0));
    
    for(int i = 0; i < n_segments; i++) {
        for(int j = i+1; j < n_segments; j++) {  
            double writhe = Gauss_Int_4_segment(
                coords[i], coords[i+1],    // First segment
                coords[j], coords[j+1]     // Second segment
            );
            writhe_matrix[i][j] = writhe;
            writhe_matrix[j][i] = writhe; 
        }
    }
    
    return writhe_matrix;
}

// Modified calculateWrithePenalty function to use the new method
double calculateWrithePenalty(ktlMolecule &mol) {
    double writhePenalty = 0.0;
    
    // Loop through each chain in the molecule
    for(int chainNum = 1; chainNum <= mol.noChains(); chainNum++) {
        // Get coordinates for this chain
        std::vector<std::vector<point>> coords = mol.getSubsecCoordinates(chainNum);
        auto flat_coords = flatten_coords(coords);
        
        // Calculate writhe matrix for this chain
        auto writhe_matrix = calculateWritheMatrix(flat_coords);
        
        // Sum up all writhe contributions
        double chainWrithe = 0.0;
        for(size_t i = 0; i < writhe_matrix.size(); i++) {
            for(size_t j = 0; j < writhe_matrix[i].size(); j++) {
                chainWrithe += writhe_matrix[i][j];
            }
        }
        
        // Get the length of this subsection
        double secLen = double(mol.getSubsecSize(chainNum));
        
        // Calculate lower bound - same formula as in original implementation
        double lowerBound = std::pow((secLen/7.5), 1.6) - 3.0;
        
        // Calculate penalty using sigmoid function
        writhePenalty += chainWrithe - lowerBound;
    }
    
    return writhePenalty;
}

// Function to compute local Average Crossing Number (lACN)
std::vector<double> compute_local_acn(const std::vector<std::vector<double>>& sigma_matrix, int window_size) {
    int N = sigma_matrix.size();
    std::vector<double> local_acn_values;
    
    // Glide window over the protein chain
    for(int start = 0; start <= N - window_size; start++) {
        int end = start + window_size - 1;
        double total = 0.0;
        
        // Sum up writhe contributions within the window
        for(int i = start; i <= end; i++) {
            for(int j = start; j < i - 1; j++) {
                total += std::abs(sigma_matrix[i][j]);
            }
        }
        
        local_acn_values.push_back(total);
    }
    
    return local_acn_values;
}


/* Updating the structure */

void refine_structure( ktlMolecule& mol, std::vector<char>& sequence, std::vector<SectionInfo>& flexible_indices, const ExperimentalData& exp_data, int noSteps ) {

    auto coords = mol.getCoordinates();
    auto flat_coords = flatten_coords(coords);

    auto CA_distances = calculate_distances(flat_coords);
    auto min_dist = find_minimum_distance(CA_distances);
    if (min_dist < 3) {  std::cout << "*Warning* input CA distance less than 3.0 A, found: " << min_dist << "A" << std::endl; }

    auto centers = process_structure(flat_coords, sequence);
    auto I_mod = calculate_saxs_implicit(centers);
    auto q_mod = init_q_values();

    auto I_mod_interp = calculate_intensity_at_experimental_q(exp_data, q_mod, I_mod);
    double scale_factor = calculate_scale_factor(exp_data, I_mod_interp);
    double current_chi_squared = calculate_chi_squared(exp_data, I_mod_interp, scale_factor);
    double current_writhe = calculateWrithePenalty(mol);
    std::cout << "Initial Chi squared: " << current_chi_squared << " (Writhe: " << current_writhe << ")" << std::endl;

    std::cout << "\nStarting refinement..." << std::endl;
    for (int i = 0; i < noSteps; i++) {
        // Update progress every 1% or when improvement found
        if (i % (noSteps/100 + 1) == 0) {
            updateProgressBar(i, noSteps);
        }

        ktlMolecule new_mol = mol;
        int rand_idx = rand() % flexible_indices.size();
        SectionInfo section_info = flexible_indices[rand_idx];
        
        new_mol.changeMoleculeSingleMulti(section_info.section_index, section_info.chain_number);
        bool violated = new_mol.checkCalphas(section_info.chain_number, mol);
        if (violated) {
            // std::cout << "Check calphas failed - skipping" << std::endl;
            continue;
        }

        auto new_coords = new_mol.getCoordinates();
        auto new_flat_coords = flatten_coords(new_coords);
        auto new_CA_distances = calculate_distances(new_flat_coords);
        auto new_min_dist = find_minimum_distance(new_CA_distances);

        if (new_min_dist < 3) {
            // std::cout << "Min CA distance less than 3.0 A - skipping" << std::endl;
            continue;
        };

        auto new_centers = process_structure(new_flat_coords, sequence);
        auto new_I_mod = calculate_saxs_implicit(new_centers);
        auto new_I_mod_interp = calculate_intensity_at_experimental_q(exp_data, q_mod, new_I_mod);
        double new_scale_factor = calculate_scale_factor(exp_data, new_I_mod_interp);
        double new_chi_squared = calculate_chi_squared(exp_data, new_I_mod_interp, new_scale_factor);
        double new_writhe = calculateWrithePenalty(new_mol);

        // Add writhe penalty to chi-squared (can adjust weight)
        double writhe_weight = 0.1;  // Adjust this weight as needed
        double writhe_penalty = writhe_weight * std::abs(new_writhe - current_writhe);
        new_chi_squared += writhe_penalty;

        if (new_chi_squared < current_chi_squared) {
            mol = new_mol;
            current_chi_squared = new_chi_squared;
            current_writhe = new_writhe;
            std::cout << "\nImproved Chi squared: " << new_chi_squared 
                     << " (Writhe: " << new_writhe << ")" << std::endl;
            updateProgressBar(i, noSteps);
        } else {
            // std::cout << "No improvement: " << new_chi_squared << std::endl;
        }
    }
    updateProgressBar(noSteps, noSteps);  // Ensure 100% at end
    std::cout << "\nRefinement complete." << std::endl;
    std::cout << "Final Chi squared: " << current_chi_squared << std::endl;
}

std::vector<int> readFlexible(const std::string& path, int structureIndex) {
    std::vector<int> flexible_sections;
    std::string filename = path + "varyingSectionSecondary" + std::to_string(structureIndex + 1) + ".dat";
    
    std::ifstream file(filename);
    if (!file.is_open()) {
        std::cerr << "Failed to open varying section file: " << filename << std::endl;
        return flexible_sections;
    }

    std::string line;
    while (std::getline(file, line) && !line.empty()) {
        std::stringstream ss(line);
        int index;
        while (ss >> index) {
            flexible_sections.push_back(index);
        }
    }
    
    file.close();
    return flexible_sections;
}

void writeCoordinatesToFile( ktlMolecule& mol, const std::string& filename ) {
    auto coords = mol.getCoordinates();
    std::ofstream outfile(filename);
    
    if (!outfile.is_open()) {
        std::cerr << "Failed to open output file: " << filename << std::endl;
        return;
    }
    
    // For each chain
    for ( auto& chain : coords) {
        // Write coordinates
        for ( auto& point : chain) {
            outfile << point.getX() << " " << point.getY() << " " << point.getZ() << std::endl;
        }
    }
    
    outfile.close();
}


int main(int argc, char* argv[]) {

    if (argc < 4) {  
        std::cerr << "Usage: " << argv[0] 
                  << " <number_of_structures> <path> <experimental_data_file>" 
                  << std::endl;
        return 1;
    }
    
    int noStructures = std::stoi(argv[1]);
    std::string path = argv[2];
    std::string data_file = argv[3];

    // Read sequence - only for first structure atm
    std::string sequenceLoc = path + "fingerPrint1.dat";
    auto sequence = extract_sequence(sequenceLoc.c_str());
    
    // Initialize and read molecular structures
    std::vector<ktlMolecule> moleculeStructures;
    readInStructures(moleculeStructures, noStructures, path);
    ktlMolecule mol = moleculeStructures[0];

    // Read flexible sections for first structure
    auto flexible_sections = readFlexible(path, 0);
    auto flexible_indices = get_flexible_indices(mol, flexible_sections);

    // Read experimental data
    auto exp_data = read_experimental_data(data_file.c_str());

    skmt s;
    std::vector<point> skmt_coords = s.getSKMTCurve(mol);
    auto writhe_matrix = calculateWritheMatrix(skmt_coords);
    
    // Write SKMT coordinates to file
    std::ofstream skmt_file(path + "skmt_coords.txt");
    if (!skmt_file.is_open()) {
        std::cerr << "Failed to open SKMT coordinates file" << std::endl;
    } else {
        for ( auto& p : skmt_coords) {
            skmt_file << p.getX() << " " << p.getY() << " " << p.getZ() << std::endl;
        }
        skmt_file.close();
    }

    // Write writhe matrix to file 
    std::ofstream writhe_file(path + "writhe_matrix.txt");
    if (!writhe_file.is_open()) {
        std::cerr << "Failed to open writhe matrix file" << std::endl;
    } else {
        for ( auto& row : writhe_matrix) {
            for ( auto& val : row) {
                writhe_file << val << " ";
            }
            writhe_file << std::endl;
        }
        writhe_file.close();
    }
 


    // Refine structure
    int noSteps = 10000;
    // refine_structure(mol, sequence, flexible_indices, exp_data, noSteps);

    // Write final coordinates
    // std::string output_file = path + "refined_coordinates.txt";
    // writeCoordinatesToFile(mol, output_file);


    
    return 0;
}
