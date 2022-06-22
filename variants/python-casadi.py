# conda create --name casadi python=3.8 -y
# conda activate casadi
# pip install casadi

import casadi
import json

def solve_opf(filename):
    x, x0, lbx, ubx, cons, lbg, ubg = [], [], [], [], [], [], []
    with open(filename, 'r') as io:
        ref = json.load(io)
    va, vm = {}, {}
    for k in ref['bus']:
        va[k] = casadi.SX.sym('va%s' % k)
        x.append(va[k])
        x0.append(0.0)
        lbx.append(-casadi.inf)
        ubx.append(casadi.inf)
        vm[k] = casadi.SX.sym('vm%s' % k)
        x.append(vm[k])
        x0.append(1.0)
        lbx.append(ref['bus'][k]['vmin'])
        ubx.append(ref['bus'][k]['vmax'])
    pg, qg = {}, {}
    for k in ref['gen']:
        pg[k] = casadi.SX.sym('pg%s' % k)
        x.append(pg[k])
        x0.append(0.0)
        lbx.append(ref['gen'][k]['pmin'])
        ubx.append(ref['gen'][k]['pmax'])
        qg[k] = casadi.SX.sym('qg%s' % k)
        x.append(qg[k])
        x0.append(0.0)
        lbx.append(ref['gen'][k]['qmin'])
        ubx.append(ref['gen'][k]['qmax'])
    p, q = {}, {}
    for arc in ref['arcs']:
        k = str(arc)
        a = ref['branch'][str(arc[0])]['rate_a']
        p[k] = casadi.SX.sym('p%s' % k)
        x.append(p[k])
        x0.append(0.0)
        lbx.append(-a)
        ubx.append(a)
        q[k] = casadi.SX.sym('q%s' % k)
        x.append(q[k])
        x0.append(0.0)
        lbx.append(-a)
        ubx.append(a)
    f = sum(
            ref['gen'][k]['cost'][0] * pg[k]**2 + 
            ref['gen'][k]['cost'][1] * pg[k] + 
            ref['gen'][k]['cost'][2] for k in ref['gen']
        )
    for k in ref['ref_buses']:
        cons.append(va[k])
        lbg.append(0)
        ubg.append(0)
    for i in ref['bus']:
        bus_loads = [ref['load'][str(l)] for l in ref['bus_loads'][i]]
        bus_shunts = [ref['shunt'][str(s)] for s in ref['bus_shunts'][i]]
        cons.append(
            sum(p[str(k)] for k in ref['bus_arcs'][i]) - 
            sum(pg[str(g)] for g in ref['bus_gens'][i]) +
            sum(load['pd'] for load in bus_loads) +
            sum(shunt['gs'] for shunt in bus_shunts) * vm[i]**2
        )
        lbg.append(0)
        ubg.append(0)
        cons.append(
            sum(q[str(k)] for k in ref['bus_arcs'][i]) - 
            sum(qg[str(g)] for g in ref['bus_gens'][i]) +
            sum(load['qd'] for load in bus_loads) -
            sum(shunt['bs'] for shunt in bus_shunts) * vm[i]**2
        )
        lbg.append(0)
        ubg.append(0)
    for i in ref['branch']:
        branch = ref['branch'][i]
        f_idx = str([int(i), branch["f_bus"], branch["t_bus"]])
        t_idx = str([int(i), branch["t_bus"], branch["f_bus"]])
        p_fr = p[f_idx]
        q_fr = q[f_idx]
        p_to = p[t_idx]
        q_to = q[t_idx]
        vm_fr = vm[str(branch["f_bus"])]
        vm_to = vm[str(branch["t_bus"])]
        va_fr = va[str(branch["f_bus"])]
        va_to = va[str(branch["t_bus"])]
        g, b = branch['g'], branch['b']
        tr, ti = branch['tr'], branch['ti']
        ttm = tr**2 + ti**2
        g_fr = branch["g_fr"]
        b_fr = branch["b_fr"]
        g_to = branch["g_to"]
        b_to = branch["b_to"]
        cons.append(
            (g+g_fr)/ttm*vm_fr**2 + (-g*tr+b*ti)/ttm*(vm_fr*vm_to*casadi.cos(va_fr-va_to)) + (-b*tr-g*ti)/ttm*(vm_fr*vm_to*casadi.sin(va_fr-va_to)) - p_fr
        )
        cons.append(
            -(b+b_fr)/ttm*vm_fr**2 - (-b*tr-g*ti)/ttm*(vm_fr*vm_to*casadi.cos(va_fr-va_to)) + (-g*tr+b*ti)/ttm*(vm_fr*vm_to*casadi.sin(va_fr-va_to)) - q_fr
        )
        cons.append(
            (g+g_to)*vm_to**2 + (-g*tr-b*ti)/ttm*(vm_to*vm_fr*casadi.cos(va_to-va_fr)) + (-b*tr+g*ti)/ttm*(vm_to*vm_fr*casadi.sin(va_to-va_fr)) - p_to
        )
        cons.append(
            -(b+b_to)*vm_to**2 - (-b*tr+g*ti)/ttm*(vm_to*vm_fr*casadi.cos(va_to-va_fr)) + (-g*tr-b*ti)/ttm*(vm_to*vm_fr*casadi.sin(va_to-va_fr)) - q_to
        )
        for i in range(4):
            lbg.append(0)
            ubg.append(0)
        cons.append(va_fr - va_to)
        lbg.append(branch['angmin'])
        ubg.append(branch['angmax'])
        cons.append(p_fr**2 + q_fr**2)
        lbg.append(-casadi.inf)
        ubg.append(branch['rate_a']**2)
        cons.append(p_to**2 + q_to**2)
        lbg.append(-casadi.inf)
        ubg.append(branch['rate_a']**2)
    model = casadi.nlpsol(
        'model', 
        'ipopt', 
        {'x': casadi.vcat(x), 'f': f, 'g': casadi.vcat(cons)},
    )
    solution = model(
        lbx = lbx,
        ubx = ubx,
        lbg = lbg,
        ubg = ubg,
        x0 = x0,
    )
    return solution

if __name__ == "__main__":
    ret = solve_opf("../data/pglib_opf_case5_pjm.json")
    print(ret)
