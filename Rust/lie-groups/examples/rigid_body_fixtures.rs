use lie_groups::{
    rigid_body::{Model, State, StepOptions, Wrench},
    Matrix,
};
fn state(s: State) -> String {
    format!(
        "{{\"r\":{:?},\"v\":{:?},\"q\":{:?},\"omega\":{:?}}}",
        s.position, s.velocity, s.quaternion, s.omega_body
    )
}
fn main() -> Result<(), Box<dyn std::error::Error>> {
    let names = [
        "free_sphere",
        "free_asymmetric",
        "free_offdiag",
        "principal",
        "forced_translation",
        "forced_torque",
        "simultaneous",
        "nonidentity",
    ];
    for (k, name) in names.iter().enumerate() {
        let mut j = [[2., 0., 0.], [0., 2., 0.], [0., 0., 2.]];
        if k == 1 || k == 3 {
            j = [[1., 0., 0.], [0., 2., 0.], [0., 0., 2.5]];
        }
        if k == 2 || k == 6 || k == 7 {
            j = [[2., 0.2, 0.1], [0.2, 3., 0.3], [0.1, 0.3, 4.]];
        }
        let mut s = State {
            position: [0.1, -0.2, 0.3],
            velocity: [0.3, 0.1, -0.2],
            omega_body: [0.2, -0.4, 0.7],
            ..State::default()
        };
        if k == 1 {
            s.omega_body = [0.2, 1., 0.3];
        }
        if k == 2 || k == 6 || k == 7 {
            s.omega_body = [0.3, -0.4, 0.8];
        }
        if k == 3 {
            s.omega_body = [0., 0., 0.7];
        }
        if k == 4 || k == 5 {
            s.omega_body = [0.; 3];
        }
        if k == 7 {
            s.quaternion = [0.3_f64.cos(), 0.3_f64.sin(), 0., 0.];
        }
        let mut w = Wrench::default();
        if k == 4 {
            w.force_body = [2., -1., 0.5];
        }
        if k == 5 {
            w.torque_body = [0.1, -0.2, 0.3];
        }
        if k == 6 {
            w.force_body = [0.2, -0.1, 0.3];
            w.torque_body = [0.03, -0.02, 0.01];
        }
        let dt = if k == 1 {
            0.005
        } else if k == 6 {
            0.002
        } else {
            0.01
        };
        let count = if k == 1 { 400 } else { 100 };
        let initial = s;
        let m = Model::new(2., Matrix::from_rows(j))?;
        for _ in 0..count {
            s = m.step(s, dt, w, StepOptions::default())?.state;
        }
        println!("{{\"name\":\"{name}\",\"mass\":2,\"inertia\":{j:?},\"dt\":{dt:?},\"steps\":{count},\"initial\":{},\"wrench\":{{\"force\":{:?},\"torque\":{:?},\"force_inertial\":{:?}}},\"final\":{},\"energy\":{:?},\"H_I\":{:?}}}",state(initial),w.force_body,w.torque_body,w.force_inertial,state(s),m.kinetic_energy(s)?,m.momentum_inertial(s)?);
    }
    Ok(())
}
