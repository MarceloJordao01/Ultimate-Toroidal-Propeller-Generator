// toroidal_propeller.scad

eps = 1/128;
$fn = 100;

// ------------------------------------------------------------
// Catmull–Rom 3D
// ------------------------------------------------------------
function catmull_rom3d(u, P0, P1, P2, P3) =
    let(u2 = u*u, u3 = u2*u)
    [
        0.5 * (2*P1[0] + (-P0[0] + P2[0]) * u +
               (2*P0[0] - 5*P1[0] + 4*P2[0] - P3[0]) * u2 +
               (-P0[0] + 3*P1[0] - 3*P2[0] + P3[0]) * u3),

        0.5 * (2*P1[1] + (-P0[1] + P2[1]) * u +
               (2*P0[1] - 5*P1[1] + 4*P2[1] - P3[1]) * u2 +
               (-P0[1] + 3*P1[1] - 3*P2[1] + P3[1]) * u3),

        0.5 * (2*P1[2] + (-P0[2] + P2[2]) * u +
               (2*P0[2] - 5*P1[2] + 4*P2[2] - P3[2]) * u2 +
               (-P0[2] + 3*P1[2] - 3*P2[2] + P3[2]) * u3)
    ];


// ------------------------------------------------------------
// Spline 3D do toroide – passa por 5 pontos
// ------------------------------------------------------------
function toroidal_path_spline_3d(
    t,
    hub_d,
    hub_height,
    blade_length,
    leading_edge_blade_width_pct,
    trailing_edge_blade_width_pct,
    leading_edge_blade_xoffset,
    trailing_edge_blade_xoffset
) =
    let(
        // offset em Y convertido de %
        leadW = blade_length * (leading_edge_blade_width_pct / 100),
        trailW = blade_length * (trailing_edge_blade_width_pct / 100),

        // raio do hub
        r  = hub_d * cos(30) / 2,

        // 5 pontos principais
        A = [ r*cos(60),  r*sin(60),   hub_height/2 ],

        B = [ blade_length * leading_edge_blade_xoffset,
              +leadW,
              hub_height/2 ],

        M = [ blade_length,
              0,
              hub_height/2 ],

        C = [ blade_length * trailing_edge_blade_xoffset,
              -trailW,
              hub_height/2 ],

        D = [ r*cos(60),  r*sin(-60),  hub_height/2 ],

        // segmentação
        segw = 1/4,
        seg =
            (t < segw)       ? 0 :
            (t < 2*segw)     ? 1 :
            (t < 3*segw)     ? 2 : 3,

        u = (seg == 0) ?  t/segw :
            (seg == 1) ? (t -     segw)/segw :
            (seg == 2) ? (t - 2 * segw)/segw :
                         (t - 3 * segw)/segw,

        P = [A, B, M, C, D],

        i0 = (seg == 0) ? 0 : seg-1,
        i1 = seg,
        i2 = seg+1,
        i3 = (seg+2 > 4) ? 4 : seg+2,

        P0 = P[i0],
        P1 = P[i1],
        P2 = P[i2],
        P3 = P[i3]
    )
    catmull_rom3d(u, P0, P1, P2, P3);


// gera lista de pontos
function toroidal_path_points(
    steps,
    hub_d,
    hub_height,
    blade_length,
    leading_edge_blade_width_pct,
    trailing_edge_blade_width_pct,
    leading_edge_blade_xoffset,
    trailing_edge_blade_xoffset
) =
    [
        for (i = [0 : steps])
            let(t = i/steps)
            toroidal_path_spline_3d(
                t, hub_d, hub_height, blade_length,
                leading_edge_blade_width_pct,
                trailing_edge_blade_width_pct,
                leading_edge_blade_xoffset,
                trailing_edge_blade_xoffset
            )
    ];


// ------------------------------------------------------------
// DEBUG leve
// ------------------------------------------------------------
module debug_toroidal_path(
    hub_d,
    hub_height,
    blade_length,
    leading_edge_blade_width_pct,
    trailing_edge_blade_width_pct,
    leading_edge_blade_xoffset,
    trailing_edge_blade_xoffset,
    steps = 32,
    radius = 0.6
){
    pts = toroidal_path_points(
        steps, hub_d, hub_height, blade_length,
        leading_edge_blade_width_pct,
        trailing_edge_blade_width_pct,
        leading_edge_blade_xoffset,
        trailing_edge_blade_xoffset
    );

    for (p = pts)
        translate(p) sphere(r = radius);
}


// ------------------------------------------------------------
// Módulo principal
// ------------------------------------------------------------
module toroidal_propeller(
    blades = 3,
    height = 6,
    blade_length = 68,
    blade_width = 42,
    blade_thickness = 4,
    blade_hole_offset = 1.4,
    blade_attack_angle = 35,
    blade_offset = -6,
    blade_safe_direction = "PREV",
    hub_height = 6,
    hub_d = 16,
    hub_screw_d = 5.5,
    hub_notch_height = 0,
    hub_notch_d = 0,

    // novos parâmetros (percentuais e offsets)
    leading_edge_blade_width = 18,       // agora interpretado como %
    trailing_edge_blade_width = 18,      // agora interpretado como %
    leading_edge_blade_xoffset = 0.25,   // fração da blade_length
    trailing_edge_blade_xoffset = 0.25
){
    difference() {
        union() {
            for (a = [0 : blades - 1])
                rotate([0,0,a*(360/blades)])
                    debug_toroidal_path(
                        hub_d, hub_height, blade_length,
                        leading_edge_blade_width,
                        trailing_edge_blade_width,
                        leading_edge_blade_xoffset,
                        trailing_edge_blade_xoffset
                    );

            rotate([0,0,30])
                cylinder(d = hub_d, h = hub_height, $fn = 6);
        }

        translate([0,0,-eps])
            cylinder(d = hub_screw_d, h = hub_height + 2*eps);

        translate([0,0,-eps])
            cylinder(d = hub_notch_d, h = hub_notch_height + eps);
    }
}


toroidal_propeller();
