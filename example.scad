use <src/toroidal_propeller.scad>
$fn = 100;                      // how polligonall you want the model

toroidal_propeller(
    blades = 1,                     // number of blades
    height = 6,                     // height
    blade_length = 68,              // blade length in mm
    blade_width = 42,               // (reservado para geometria futura)
    blade_thickness = 4,            // (reservado para geometria futura)
    blade_hole_offset = 1.4,        // (reservado para geometria futura)
    blade_attack_angle = 35,        // (reservado para geometria futura)
    blade_offset = -6,              // (reservado para geometria futura)
    blade_safe_direction = "PREV",  // (reservado para geometria futura)
    hub_height = 6,                 // Hub height
    hub_d = 16,                     // hub diameter (círculo que circunscreve o hexágono)
    hub_screw_d = 5.5,              // hub screw diameter
    hub_notch_height = 0,           // height for the notch 
    hub_notch_d = 0,                // diameter for the notch
    // --- NOVOS PARÂMETROS DE FORMA DO CAMINHO ---
    leading_edge_blade_width = 18,       // agora interpretado como %
    trailing_edge_blade_width = 18,      // agora interpretado como %
    leading_edge_blade_xoffset = 50,   // fração da blade_length
    trailing_edge_blade_xoffset = 25,  // fração da blade_length
    // perfis NACA:
    naca_profile_1 = [2,4,1,2],  // NACA 2412
    naca_profile_2 = [0,0,1,2],  // NACA 0012
    use_profile    = 2,           // alterna entre 1 e 2 pra ver a diferença
    path_portion = 1.0
);
