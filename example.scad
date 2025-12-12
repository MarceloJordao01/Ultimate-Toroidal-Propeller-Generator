use <src/toroidal_propeller.scad>
$fn = 100;                      // how polligonall you want the model

toroidal_propeller(
    blades = 2,                         // number of blades
    hub_height = 6,                     // Hub height
    hub_d = 16,                         // hub diameter (círculo que circunscreve o hexágono)
    hub_screw_d = 5.5,                  // hub screw diameter
    hub_notch_height = 0,               // height for the notch 
    hub_notch_d = 0,                    // diameter for the notch
    // --- NOVOS PARÂMETROS DE FORMA DO CAMINHO ---
    blade_offset = 2,                   // defazagem das laminas em Z
    leading_edge_blade_width = 18,      // distancia R1 do path do toroide, porcentagem do valor de blade_length/2
    trailing_edge_blade_width = 18,     // distancia R2 do path do toroide, porcentagem do valor de blade_length/2
    leading_edge_blade_xoffset = 50,    // distancia R1x do path do toroide, porcentagem do valor de blade_length/2
    trailing_edge_blade_xoffset = 60,   // distancia R2x do path do toroide, porcentagem do valor de blade_length/2
    // perfis NACA:
    naca_profiles=["2412","2412"],      // perfil inicial e final por enquanto
    chords=[6, 6],                      // tamanho das chords
    attack_angles=[15, -10],            // angulo de ataque
    chord_pivot_pcts = [-10, 0],        // posicao do 
    path_portion = 1.0
);
