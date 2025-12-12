// toroidal_propeller.scad
include <BOSL2/std.scad>;

eps = 1/128;
$fn = 100;

// ============================================================
// Helpers
// ============================================================

// "2412" -> [2,4,1,2]
function naca_digits_from_string(s) =
    [ for (i=[0:len(s)-1]) ord(s[i]) - ord("0") ];

function list_get_or_last(L, i, fallback=undef) =
    (L == undef || len(L)==0) ? fallback :
    (i < len(L)) ? L[i] : L[len(L)-1];

function lerp(a,b,u) = a + (b-a)*u;
function lerp2(p,q,u) = [ lerp(p[0],q[0],u), lerp(p[1],q[1],u) ];

// 2D rotation
function rot2d_pt(p, ang) =
    let(c=cos(ang), s=sin(ang))
    [ p[0]*c - p[1]*s, p[0]*s + p[1]*c ];

function rot2d_path(path2d, ang) =
    [ for (p=path2d) rot2d_pt(p, ang) ];

// 2D -> 3D in XY plane
function path2d_to3d(path2d) =
    [ for (p=path2d) [p[0], p[1], 0] ];

// Safe normalize
function vnorm(v) = norm(v);
function vunit(v) = (vnorm(v) < 1e-9) ? [0,0,1] : (v / vnorm(v));

// Cross
function vcross(a,b) = cross(a,b);

// Rotate vector v around axis t by angle deg (Rodrigues). t must be unit.
function rot_about_axis(v, t, ang) =
    let(c = cos(ang), s = sin(ang))
    v*c + vcross(t, v)*s + t*(t*v)*(1-c);   // (t*v) is dot(t,v)

// Build 4x4 transform matrix from basis vectors and origin.
// Local axes: X=n, Y=b, Z=t, origin=p
function frame_matrix(p, n, b, t) =
    [
        [ n[0], b[0], t[0], p[0] ],
        [ n[1], b[1], t[1], p[1] ],
        [ n[2], b[2], t[2], p[2] ],
        [ 0,    0,    0,    1    ]
    ];

// Tangent at index i from path list
function path_tangent(path, i) =
    let(
        N = len(path),
        p0 = (i==0)   ? path[0]   : path[i-1],
        p1 = (i==N-1) ? path[N-1] : path[i+1]
    )
    vunit(p1 - p0);

// Make a stable normal/binormal given tangent and a reference up
function make_nb_from_tangent(t, up=[0,0,1]) =
    let(
        up2 = (abs(t*up) > 0.95) ? [1,0,0] : up,
        b = vunit(vcross(t, up2)),
        n = vunit(vcross(b, t))
    )
    [n, b];


// ============================================================
// Catmull–Rom 3D
// ============================================================
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


// ============================================================
// Spline 3D do path toroidal (A–B–M–C–D)
// derivadas controladas em A e D via pontos fantasmas
// ============================================================
function toroidal_path_spline_3d(
    t,
    hub_d,
    hub_height,
    blade_length,
    blade_offset,
    leading_edge_blade_width_pct,
    trailing_edge_blade_width_pct,
    leading_edge_blade_xoffset_pct,
    trailing_edge_blade_xoffset_pct,
    angle_A = 60,
    angle_D = -240
) =
    let(
        leadW  = blade_length * (leading_edge_blade_width_pct  / 100),
        trailW = blade_length * (trailing_edge_blade_width_pct / 100),

        leadX  = blade_length * (leading_edge_blade_xoffset_pct  / 100),
        trailX = blade_length * (trailing_edge_blade_xoffset_pct / 100),

        r = hub_d * cos(30) / 2,

        A = [ r*cos(60),  r*sin(60),    hub_height/2 + blade_offset/2 ],
        B = [ leadX,      +leadW,       hub_height/2 ],
        M = [ blade_length, 0,          hub_height/2 ],
        C = [ trailX,     -trailW,      hub_height/2 ],
        D = [ r*cos(60),  r*sin(-60),   hub_height/2 - blade_offset/2 ],

        scale_A = norm(B - A),
        scale_D = norm(D - C),

        Ta = [ scale_A*cos(angle_A), scale_A*sin(angle_A), 0 ],
        Td = [ scale_D*cos(angle_D), scale_D*sin(angle_D), 0 ],

        A0 = B - 2*Ta,
        D3 = C + 2*Td,

        segw = 1/4,
        seg =
            (t < segw)   ? 0 :
            (t < 2*segw) ? 1 :
            (t < 3*segw) ? 2 : 3,

        u = (seg == 0) ?  t/segw :
            (seg == 1) ? (t - segw)/segw :
            (seg == 2) ? (t - 2*segw)/segw :
                         (t - 3*segw)/segw,

        P0 = (seg == 0) ? A0 :
             (seg == 1) ? A  :
             (seg == 2) ? B  : M,

        P1 = (seg == 0) ? A  :
             (seg == 1) ? B  :
             (seg == 2) ? M  : C,

        P2 = (seg == 0) ? B  :
             (seg == 1) ? M  :
             (seg == 2) ? C  : D,

        P3 = (seg == 0) ? M  :
             (seg == 1) ? C  :
             (seg == 2) ? D  : D3
    )
    catmull_rom3d(u, P0, P1, P2, P3);


// ============================================================
// Lista de pontos do path toroidal
// ============================================================
function toroidal_path_points(
    steps,
    hub_d,
    hub_height,
    blade_length,
    blade_offset,
    leading_edge_blade_width_pct,
    trailing_edge_blade_width_pct,
    leading_edge_blade_xoffset_pct,
    trailing_edge_blade_xoffset_pct,
    t_start = 0,
    t_end   = 1
) =
[
    for (i = [0:steps])
        let(
            u = i/steps,
            t = t_start + (t_end - t_start) * u
        )
        toroidal_path_spline_3d(
            t,
            hub_d,
            hub_height,
            blade_length,
            blade_offset,
            leading_edge_blade_width_pct,
            trailing_edge_blade_width_pct,
            leading_edge_blade_xoffset_pct,
            trailing_edge_blade_xoffset_pct
        )
];


// ============================================================
// Perfil NACA 4 dígitos (2D) - fechado
// ============================================================
function naca4_path(digits, chord=8, pts=40) =
    let(
        m = digits[0]/100,
        p = digits[1]/10,
        tt = (digits[2]*10 + digits[3]) / 100,

        upper = [
            for (i=[0:pts])
                let(
                    xc = i/pts,
                    yt = 5*tt*(0.2969*sqrt(xc)-0.1260*xc-0.3516*xc*xc
                              +0.2843*xc*xc*xc-0.1015*xc*xc*xc*xc),
                    yc = (m==0)?0:(xc<p ? m*xc*xc/(p*p)*(2*p-xc)
                                       : m*(1-xc)*(1-xc)/((1-p)*(1-p))*(1+2*p-xc-2*p)),
                    dyc = (m==0)?0:(xc<p ? 2*m/(p*p)*(p-xc)
                                         : 2*m/((1-p)*(1-p))*(p-xc)),
                    th = atan(dyc)
                )
                [(xc-yt*sin(th))*chord, (yc+yt*cos(th))*chord]
        ],

        lower = [
            for (i=[pts:-1:0])
                let(
                    xc = i/pts,
                    yt = 5*tt*(0.2969*sqrt(xc)-0.1260*xc-0.3516*xc*xc
                              +0.2843*xc*xc*xc-0.1015*xc*xc*xc*xc),
                    yc = (m==0)?0:(xc<p ? m*xc*xc/(p*p)*(2*p-xc)
                                       : m*(1-xc)*(1-xc)/((1-p)*(1-p))*(1+2*p-xc-2*p)),
                    dyc = (m==0)?0:(xc<p ? 2*m/(p*p)*(p-xc)
                                         : 2*m/((1-p)*(1-p))*(p-xc)),
                    th = atan(dyc)
                )
                [(xc+yt*sin(th))*chord, (yc-yt*cos(th))*chord]
        ],

        closed = concat(upper, lower)
    )
    [
        for (pp = closed)
            [(pp[0]-0.25*chord), pp[1]]
    ];


// ============================================================
// Asa NACA ao longo do path usando SKIN
// ============================================================
module toroidal_wing_naca_skin(
    hub_d,
    hub_height,
    blade_length,
    blade_offset,
    leading_edge_blade_width_pct,
    trailing_edge_blade_width_pct,
    leading_edge_blade_xoffset_pct,
    trailing_edge_blade_xoffset_pct,

    naca_profiles,     // usa 0 e 1
    chords,            // usa 0 e 1
    attack_angles,     // usa 0 e 1
    chord_pivot_pcts,

    naca_pts,
    path_steps,
    path_portion,
    extra_twist_deg
){
    if (naca_profiles == undef || len(naca_profiles) < 2) {
        echo("ERRO: naca_profiles precisa ter pelo menos 2 termos.");
    } else {

        path = toroidal_path_points(
            path_steps,
            hub_d,
            hub_height,
            blade_length,
            blade_offset,
            leading_edge_blade_width_pct,
            trailing_edge_blade_width_pct,
            leading_edge_blade_xoffset_pct,
            trailing_edge_blade_xoffset_pct,
            0,
            path_portion
        );

        digits0 = naca_digits_from_string(naca_profiles[0]);
        digits1 = naca_digits_from_string(naca_profiles[1]);

        chord0 = list_get_or_last(chords, 0, 4);
        chord1 = list_get_or_last(chords, 1, chord0);

        atk0   = list_get_or_last(attack_angles, 0, 0);
        atk1   = list_get_or_last(attack_angles, 1, atk0);

        pivot0 = list_get_or_last(chord_pivot_pcts, 0, 25) / 100;
        pivot1 = list_get_or_last(chord_pivot_pcts, 1, pivot0*100) / 100;   

        prof0_base = naca4_path(digits0, chord0, naca_pts);
        prof1_base = naca4_path(digits1, chord1, naca_pts);

        sections =
        [
            for (i=[0:len(path)-1])
                let(
                    u = (len(path)==1) ? 0 : i/(len(path)-1),

                    chord_u = lerp(chord0, chord1, u),
                    atk_u   = lerp(atk0,   atk1,   u),

                    s0 = (chord0==0) ? 1 : (chord_u/chord0),
                    s1 = (chord1==0) ? 1 : (chord_u/chord1),

                    p0 = [ for (p=prof0_base) [p[0]*s0, p[1]*s0] ],
                    p1 = [ for (p=prof1_base) [p[0]*s1, p[1]*s1] ],

                    prof_u = [ for (k=[0:len(p0)-1]) lerp2(p0[k], p1[k], u) ],

                    pivot_u = lerp(pivot0, pivot1, u),
                    dx = (pivot_u - 0.25) * chord_u,
                    prof_u_pivot = [ for (p=prof_u) [p[0] + dx, p[1]] ],
                    prof_u_rot = rot2d_path(prof_u_pivot, atk_u + 90),

                    t = path_tangent(path, i),
                    nb = make_nb_from_tangent(t, [0,0,1]),
                    n0 = nb[0],
                    b0 = nb[1],

                    tw = u * extra_twist_deg,
                    n = rot_about_axis(n0, t, tw),
                    b = rot_about_axis(b0, t, tw),

                    T = frame_matrix(path[i], n, b, t),

                    sec3d = apply(T, path2d_to3d(prof_u_rot))
                )
                sec3d
        ];

        // >>> CORREÇÃO PARA SUA VERSÃO DO BOSL2:
        skin(sections, slices=1, caps=true);
    }
}


// ============================================================
// Módulo principal
// ============================================================
module toroidal_propeller(
    blades = 3,
    blade_length = 68,
    blade_offset = -6,

    hub_height = 6,
    hub_d = 16,
    hub_screw_d = 5.5,
    hub_notch_height = 0,
    hub_notch_d = 0,

    leading_edge_blade_width = 18,
    trailing_edge_blade_width = 18,
    leading_edge_blade_xoffset = 25,
    trailing_edge_blade_xoffset = 25,

    naca_profiles = ["2412","0012"],
    chords = [6, 3.5],
    attack_angles = [8, 2],
    chord_pivot_pcts = [25, 25],

    naca_pts    = 80,
    path_steps  = 64,
    path_portion = 1.0,
    extra_twist_deg = -180
){
    difference() {
        union() {
            for (a=[0:blades-1])
                rotate([0,0,a*360/blades])
                    toroidal_wing_naca_skin(
                        hub_d, hub_height,
                        blade_length, blade_offset,
                        leading_edge_blade_width,
                        trailing_edge_blade_width,
                        leading_edge_blade_xoffset,
                        trailing_edge_blade_xoffset,

                        naca_profiles,
                        chords,
                        attack_angles,
                        chord_pivot_pcts,

                        naca_pts,
                        path_steps,
                        path_portion,
                        extra_twist_deg
                    );

            rotate([0,0,30])
                cylinder(d=hub_d, h=hub_height, $fn=6);
        }

        translate([0,0,-eps])
            cylinder(d=hub_screw_d, h=hub_height+2*eps);

        if (hub_notch_height > 0 && hub_notch_d > 0)
            translate([0,0,-eps])
                cylinder(d=hub_notch_d, h=hub_notch_height+eps);
    }
}

toroidal_propeller();
