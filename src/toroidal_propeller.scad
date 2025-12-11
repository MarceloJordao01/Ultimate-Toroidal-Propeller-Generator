// toroidal_propeller.scad

// BOSL2
include <BOSL2/std.scad>;

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
// leading/trailing xoffset em % (0–100)
// ------------------------------------------------------------
function toroidal_path_spline_3d(
    t,
    hub_d,
    hub_height,
    blade_length,
    leading_edge_blade_width_pct,
    trailing_edge_blade_width_pct,
    leading_edge_blade_xoffset_pct,
    trailing_edge_blade_xoffset_pct
) =
    let(
        // offsets em Y convertidos de %
        leadW  = blade_length * (leading_edge_blade_width_pct  / 100),
        trailW = blade_length * (trailing_edge_blade_width_pct / 100),

        // offsets em X convertidos de %
        leadX  = blade_length * (leading_edge_blade_xoffset_pct  / 100),
        trailX = blade_length * (trailing_edge_blade_xoffset_pct / 100),

        // raio do hub
        r  = hub_d * cos(30) / 2,

        // 5 pontos principais
        A = [ r*cos(60),  r*sin(60),   hub_height/2 ],

        B = [ leadX,
              +leadW,
              hub_height/2 ],

        M = [ blade_length,
              0,
              hub_height/2 ],

        C = [ trailX,
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


// ------------------------------------------------------------
// gera lista de pontos do path toroidal, com t_start/t_end
// ------------------------------------------------------------
function toroidal_path_points(
    steps,
    hub_d,
    hub_height,
    blade_length,
    leading_edge_blade_width_pct,
    trailing_edge_blade_width_pct,
    leading_edge_blade_xoffset_pct,
    trailing_edge_blade_xoffset_pct,
    t_start = 0,
    t_end   = 1
) =
    [
        for (i = [0 : steps])
            let(
                u = i/steps,
                t = t_start + (t_end - t_start) * u
            )
            toroidal_path_spline_3d(
                t,
                hub_d,
                hub_height,
                blade_length,
                leading_edge_blade_width_pct,
                trailing_edge_blade_width_pct,
                leading_edge_blade_xoffset_pct,
                trailing_edge_blade_xoffset_pct
            )
    ];


// ------------------------------------------------------------
// Geração de perfil NACA 4 dígitos
// digits = [d1,d2,d3,d4] -> NACA d1 d2 d3 d4  (ex: [2,4,1,2] = 2412)
// chord fixado em 8 mm para manter dimensão máxima ~8
// Perfil é espelhado em X (apontando pro lado desejado).
// ------------------------------------------------------------
function naca4_path(digits, chord=8, pts=80) =
    let(
        m = digits[0]/100,
        p = digits[1]/10,
        t = (digits[2]*10 + digits[3]) / 100,

        // upper surface (LE -> TE)
        upper = [
            for (i=[0:pts])
                let(
                    xc = i/pts,
                    x  = xc,
                    yt = 5*t*(
                          0.2969*sqrt(xc)
                        - 0.1260*xc
                        - 0.3516*xc*xc
                        + 0.2843*xc*xc*xc
                        - 0.1015*xc*xc*xc*xc
                    ),
                    yc = (m==0) ? 0 :
                         (xc < p ?
                            m*pow(xc,2)/(p*p) * (2*p - xc) :
                            m*pow(1-xc,2)/pow(1-p,2) * (1 + 2*p - xc - 2*p)
                         ),
                    dyc_dx = (m==0) ? 0 :
                             (xc < p ?
                                2*m/p/p*(p - xc) :
                                2*m/pow(1-p,2)*(p - xc)
                             ),
                    theta = atan(dyc_dx),
                    xu = x - yt*sin(theta),
                    yu = yc + yt*cos(theta)
                )
                [xu*chord, yu*chord]
        ],

        // lower surface (TE -> LE)
        lower = [
            for (i=[pts:-1:0])
                let(
                    xc = i/pts,
                    x  = xc,
                    yt = 5*t*(
                          0.2969*sqrt(xc)
                        - 0.1260*xc
                        - 0.3516*xc*xc
                        + 0.2843*xc*xc*xc
                        - 0.1015*xc*xc*xc*xc
                    ),
                    yc = (m==0) ? 0 :
                         (xc < p ?
                            m*pow(xc,2)/(p*p) * (2*p - xc) :
                            m*pow(1-xc,2)/pow(1-p,2) * (1 + 2*p - xc - 2*p)
                         ),
                    dyc_dx = (m==0) ? 0 :
                             (xc < p ?
                                2*m/p/p*(p - xc) :
                                2*m/pow(1-p,2)*(p - xc)
                             ),
                    theta = atan(dyc_dx),
                    xl = x + yt*sin(theta),
                    yl = yc - yt*cos(theta)
                )
                [xl*chord, yl*chord]
        ],

        poly_raw = concat(upper, lower),

        // move origem para 25% da corda na linha média
        poly_shifted = [
            for (p = poly_raw)
                [p[0] - 0.25*chord, p[1]]
        ],

        // espelha em X pra inverter a direção
        poly_flipped = [
            for (p = poly_shifted)
                [-p[0], p[1]]
        ]
    )
    poly_flipped;


// ------------------------------------------------------------
// Asinha NACA varrida ao longo do path toroidal (BOSL2 path_sweep)
// path_portion = 0.25 -> só 25% do caminho (debug)
// extra_twist_deg corrige a orientação da ponta
// ------------------------------------------------------------
module toroidal_wing_naca(
    hub_d,
    hub_height,
    blade_length,
    // parâmetros da spline
    leading_edge_blade_width_pct,
    trailing_edge_blade_width_pct,
    leading_edge_blade_xoffset_pct,
    trailing_edge_blade_xoffset_pct,

    // parâmetros do perfil
    naca_digits = [2,4,1,2],   // NACA 2412 por padrão
    chord       = 8,           // corda fixa em 8 mm
    naca_pts    = 80,          // resolução do perfil
    path_steps  = 64,          // resolução do caminho
    path_portion = 0.25,       // fração do caminho (0–1)
    extra_twist_deg = -180     // correção de twist pra ponta
){
    // perfil NACA em 2D (XY)
    profile = naca4_path(naca_digits, chord=chord, pts=naca_pts);

    // caminho 3D do toroide, só até path_portion
    path = toroidal_path_points(
        steps = path_steps,
        hub_d = hub_d,
        hub_height = hub_height,
        blade_length = blade_length,
        leading_edge_blade_width_pct = leading_edge_blade_width_pct,
        trailing_edge_blade_width_pct = trailing_edge_blade_width_pct,
        leading_edge_blade_xoffset_pct = leading_edge_blade_xoffset_pct,
        trailing_edge_blade_xoffset_pct = trailing_edge_blade_xoffset_pct,
        t_start = 0,
        t_end   = path_portion
    );

    // varre o perfil ao longo desse pedaço do path
    path_sweep(
        shape = profile,
        path  = path,
        caps  = true,
        twist = extra_twist_deg
    );
}


// ------------------------------------------------------------
// Módulo principal
// ------------------------------------------------------------
module toroidal_propeller(
    blades = 3,
    height = 6,
    blade_length = 68,
    blade_width = 42,      // não liga direto na corda (corda fixa em 8)
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

    // parâmetros do path (percentuais e offsets em %)
    leading_edge_blade_width = 18,      // % da blade_length em Y
    trailing_edge_blade_width = 18,     // % da blade_length em Y
    leading_edge_blade_xoffset = 25,    // % da blade_length em X
    trailing_edge_blade_xoffset = 25,   // % da blade_length em X

    // --- dois perfis NACA configuráveis ---
    naca_profile_1 = [2,4,1,2],   // NACA 2412
    naca_profile_2 = [0,0,1,2],   // NACA 0012
    use_profile    = 1,           // 1 = usa naca_profile_1, 2 = usa naca_profile_2

    naca_pts    = 80,
    path_steps  = 64,
    path_portion = 0.25,          // ainda só 25% do caminho (debug)
    extra_twist_deg = -180        // correção da ponta (inversão)
){
    chosen_profile =
        (use_profile == 1) ? naca_profile_1 :
        (use_profile == 2) ? naca_profile_2 :
                             naca_profile_1;

    difference() {
        union() {
            // asas NACA toroidais
            for (a = [0 : blades - 1])
                rotate([0,0,a*(360/blades)])
                    toroidal_wing_naca(
                        hub_d        = hub_d,
                        hub_height   = hub_height,
                        blade_length = blade_length,

                        leading_edge_blade_width_pct  = leading_edge_blade_width,
                        trailing_edge_blade_width_pct = trailing_edge_blade_width,
                        leading_edge_blade_xoffset_pct  = leading_edge_blade_xoffset,
                        trailing_edge_blade_xoffset_pct = trailing_edge_blade_xoffset,

                        naca_digits = chosen_profile,
                        chord       = 8,
                        naca_pts    = naca_pts,
                        path_steps  = path_steps,
                        path_portion = path_portion,
                        extra_twist_deg = extra_twist_deg
                    );

            // hub hexagonal
            rotate([0,0,30])
                cylinder(d = hub_d, h = hub_height, $fn = 6);
        }

        // furo central
        translate([0,0,-eps])
            cylinder(d = hub_screw_d, h = hub_height + 2*eps);

        // notch opcional
        translate([0,0,-eps])
            cylinder(d = hub_notch_d, h = hub_notch_height + eps);
    }
}


// chamada default para teste rápido
toroidal_propeller();
