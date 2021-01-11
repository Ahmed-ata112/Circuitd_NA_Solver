import cmath

# This is a sample Python script.

# omega is a global const variable
from sympy import symbols, solve

omega = 0.0
voltage_sources_count = 0  # start with 0


class node:
    # static variable
    count = 1

    def __init__(self, id=0, voltage=complex(0.0, 0.0)):
        self.id = node.count  #To calculate no of nodes
        self.components = []  #like array of string in c++ i think
        self.voltage = voltage
        self.sign = []
        node.count += 1

    def addComp(self, comp, sign):
        self.components.append(comp)
        self.sign.append(sign)


class Res:
    def __init__(self, name, res, node1, node2):  #string  int(value)  int  int
        self.name = name
        self.resistance = res
        self.first_node = node1
        self.second_node = node2
        self.admittance = complex(1.0) / res


class Cap:
    def __init__(self, name, capacitance, node1, node2):
        self.name = name
        self.capacitance = capacitance
        self.first_node = node1
        self.second_node = node2
        self.admittance = complex(0, 1) * capacitance * omega


class Ind:
    def __init__(self, name, inductance, node1, node2):
        self.name = name
        self.inductance = inductance
        self.first_node = node1
        self.second_node = node2
        self.admittance = (complex(0, 1) * inductance * omega) ** (-1)


class Vsrc:
    def __init__(self, name, mag, phase, pos_node, neg_node):
        self.name = name
        global voltage_sources_count
        self.id = voltage_sources_count
        voltage_sources_count += 1
        self.magnitude = mag
        self.phase = phase
        self.pos_terminal = pos_node
        self.neg_terminal = neg_node


class Isrc:
    def __init__(self, name, mag, phase, pos_node, neg_node):
        self.name = name
        self.magnitude = mag
        self.phase = phase
        self.pos_terminal = pos_node
        self.neg_terminal = neg_node


class Cccs:
    def __init__(self, name, coeff, passing_component, pos_node, neg_node,
                 control_pos, control_neg):
        self.name = name
        self.pos_terminal = pos_node
        self.neg_terminal = neg_node
        self.control_pos = control_pos
        self.control_neg = control_neg
        self.coefficient = coeff
        self.passing_component = passing_component


class Ccvs:
    def __init__(self, name, coeff, passing_component, pos_node, neg_node,
                 control_pos, control_neg):
        self.name = name
        global voltage_sources_count
        self.id = voltage_sources_count
        voltage_sources_count += 1
        self.pos_terminal = pos_node
        self.neg_terminal = neg_node
        self.control_pos = control_pos
        self.control_neg = control_neg
        self.coefficient = coeff
        self.passing_component = passing_component  # thats the name of the component


class Vcvs:
    def __init__(self, name, coeff, pos_node, neg_node,
                 control_pos, control_neg):
        self.name = name
        global voltage_sources_count
        self.id = voltage_sources_count
        voltage_sources_count += 1

        self.pos_terminal = pos_node
        self.neg_terminal = neg_node
        self.control_pos = control_pos
        self.control_neg = control_neg
        self.coefficient = coeff


class Vccs:
    def __init__(self, name, coeff, pos_node, neg_node,
                 control_pos, control_neg):
        self.name = name

        self.pos_terminal = pos_node
        self.neg_terminal = neg_node
        self.control_pos = control_pos
        self.control_neg = control_neg
        self.coefficient = coeff


if __name__ == '__main__':   #I dont understand this line
    resistances = []
    inductances = []
    capacitances = []

    Vsrcs = []
    Isrcs = []
    Vcvss = []
    Ccvss = []
    Vccss = []
    Cccss = []
    f = open("examples/input.txt")
    input_str = (f.read()).split() ##the split fn returns a list contains strings, when find a white space that means it is new string
## so input_str is a list of string
    no_nodes = 0
    no_voltage_sources = 0
    name_of_variables = []

    # counting the number of nodes and sources
    for i in range(len(input_str)): # loop of the number of elements inside the list to calculate no of nodes and voltage sources
        s = input_str[i]
        if s == "res" or s == "vsrc" or s == "isrc" or s == "ind" or s == "cap" or s == "vcvs" or s == "ccvs" or s == "vccs" or s == "cccs":
            pos_terminal = int(input_str[i + 2])
            neg_terminal = int(input_str[i + 3])
            no_nodes = max(no_nodes, neg_terminal, pos_terminal)

        if s == "vsrc" or s == "vcvs" or s == "ccvs":  # any voltage source
            no_voltage_sources += 1

    ########################################
    total_no_nodes = no_nodes + 1  # with the ground node
    nodes = [node() for i in range(total_no_nodes)] # creating a list of objects from class node

    V0 = symbols("V0")
    equations = [V0 for i in range(total_no_nodes)] # I think i take a default value 0 at first # i think this line initialise all equations with V0

    nodes_Voltages_sym = list(symbols("V0:%d" % (no_nodes + 1)))  # V1 V2 V3 for 3 nodes in circuit
    current_in_vs_sym = list(symbols("iV1:%d" % (no_voltage_sources + 1)))
    # equations[0] = Eq(V0, 0)
    # filling the data of the component
    i = 0
    while i < len(input_str):

        s = input_str[i]
        if s == "w":
            i += 1
            omega = float(input_str[i])

            i += 1
            continue

        if s == "res":
            i += 1
            name = input_str[i]
            i += 1
            pos_terminal = int(input_str[i])
            i += 1
            neg_terminal = int(input_str[i])
            i += 1
            val = float(input_str[i])

            r = Res(name, val, pos_terminal, neg_terminal)
            nodes[pos_terminal].addComp(r, 1)
            nodes[neg_terminal].addComp(r, -1)
            resistances.append(r)

        if s == "vsrc":
            i += 1
            name = input_str[i]
            i += 1
            pos_terminal = int(input_str[i])
            i += 1
            neg_terminal = int(input_str[i])
            i += 1
            magnitude = float(input_str[i])
            i += 1
            phase = float(input_str[i]) * cmath.pi / 180.0  # convert to radian

            v = Vsrc(name, magnitude, phase, pos_terminal, neg_terminal)
            nodes[pos_terminal].addComp(v, 1)
            nodes[neg_terminal].addComp(v, -1)
            Vsrcs.append(v)

        if s == "isrc":
            i += 1
            name = input_str[i]
            i += 1
            pos_terminal = int(input_str[i])
            i += 1
            neg_terminal = int(input_str[i])
            i += 1
            magnitude = float(input_str[i])
            i += 1
            phase = float(input_str[i]) * cmath.pi / 180.0

            isr = Isrc(name, magnitude, phase, pos_terminal, neg_terminal)
            nodes[pos_terminal].addComp(isr, 1)
            nodes[neg_terminal].addComp(isr, -1)

            Isrcs.append(isr)

        if s == "cap":
            i += 1
            name = input_str[i]
            i += 1
            pos_terminal = int(input_str[i])
            i += 1
            neg_terminal = int(input_str[i])
            i += 1
            capacitance = int(input_str[i])

            cap = Cap(name, capacitance, pos_terminal, neg_terminal)
            nodes[pos_terminal].addComp(cap, 1)
            nodes[neg_terminal].addComp(cap, -1)
            capacitances.append(cap)

        if s == "ind":
            i += 1
            name = input_str[i]
            i += 1
            pos_terminal = int(input_str[i])
            i += 1
            neg_terminal = int(input_str[i])
            i += 1
            inductance = int(input_str[i])

            ind = Ind(name, inductance, pos_terminal, neg_terminal)
            nodes[pos_terminal].addComp(ind, 1)
            nodes[neg_terminal].addComp(ind, -1)
            inductances.append(ind)

        if s == "cccs":
            i += 1
            name = input_str[i]
            i += 1
            pos_terminal = int(input_str[i])
            i += 1
            neg_terminal = int(input_str[i])
            i += 1
            control_neg = int(input_str[i])
            i += 1
            control_pos = int(input_str[i])
            i += 1
            passing_comp = input_str[i]
            i += 1
            coeff = int(input_str[i])

            cccs = Cccs(name, coeff,  passing_comp, pos_terminal,  neg_terminal, control_pos, control_neg)
            nodes[pos_terminal].addComp(cccs, 1)
            nodes[neg_terminal].addComp(cccs, -1)
            Cccss.append(cccs)

        if s == "ccvs":
            i += 1
            name = input_str[i]
            i += 1
            pos_terminal = int(input_str[i])
            i += 1
            neg_terminal = int(input_str[i])
            i += 1
            control_neg = int(input_str[i])
            i += 1
            control_pos = int(input_str[i])
            i += 1
            passing_comp = input_str[i]
            i += 1
            coeff = int(input_str[i])

            ccvs = (name, coeff, passing_comp, pos_terminal, neg_terminal, control_pos, control_neg)
            nodes[pos_terminal].addComp(ccvs, 1)
            nodes[neg_terminal].addComp(ccvs, -1)
            Ccvss.append(ccvs)

        if s == "vcvs":
            i += 1
            name = input_str[i]
            i += 1
            pos_terminal = int(input_str[i])
            i += 1
            neg_terminal = int(input_str[i])
            i += 1
            control_pos = int(input_str[i])
            i += 1
            control_neg = int(input_str[i])
            i += 1
            coeff = int(input_str[i])

            vcvs = (name, coeff, pos_terminal, neg_terminal, control_pos, control_neg)
            nodes[pos_terminal].addComp(vcvs, 1)
            nodes[neg_terminal].addComp(vcvs, -1)
            Vcvss.append(vcvs)

        if s == "vccs":
            i += 1
            name = input_str[i]
            i += 1
            pos_terminal = int(input_str[i])
            i += 1
            neg_terminal = int(input_str[i])
            i += 1
            control_pos = int(input_str[i])
            i += 1
            control_neg = int(input_str[i])
            i += 1
            coeff = int(input_str[i])

            vccs = (name, coeff, pos_terminal, neg_terminal, control_pos, control_neg)
            nodes[pos_terminal].addComp(vccs, 1)
            nodes[neg_terminal].addComp(vccs, -1)
            Vccss.append(vccs)

        i += 1

    # now every component has its values

    for r in resistances:
        V_first_node = nodes_Voltages_sym[r.first_node]  # V1
        V_sec_node = nodes_Voltages_sym[r.second_node]  # V2
        if r.first_node != 0:
            equations[r.first_node] += (V_first_node - V_sec_node) * r.admittance
        if r.second_node != 0:
            equations[r.second_node] += (V_sec_node - V_first_node) * r.admittance

    # filling the data of the I_src
    for cs in Isrcs:
        V_p_node = nodes_Voltages_sym[cs.pos_terminal]
        V_n_node = nodes_Voltages_sym[cs.neg_terminal]
        if cs.pos_terminal != 0:
            equations[cs.pos_terminal] -= cmath.rect(cs.magnitude, cs.phase)
        if cs.neg_terminal != 0:
            equations[cs.neg_terminal] += cmath.rect(cs.magnitude, cs.phase)

    for vs in Vsrcs:
        V_p_node = nodes_Voltages_sym[vs.pos_terminal]
        V_n_node = nodes_Voltages_sym[vs.neg_terminal]
        equations.append(V_p_node - V_n_node - cmath.rect(vs.magnitude, vs.phase))

        if vs.pos_terminal != 0:
            equations[vs.pos_terminal] += current_in_vs_sym[vs.id]
        if vs.neg_terminal != 0:
            equations[vs.neg_terminal] -= current_in_vs_sym[vs.id]

    for cccs in Cccss:
        V_p_node = nodes_Voltages_sym[cccs.pos_terminal]
        V_n_node = nodes_Voltages_sym[cccs.neg_terminal]
        V_control_p = nodes_Voltages_sym[cccs.control_pos]
        V_control_n = nodes_Voltages_sym[cccs.control_neg]
        if cccs.pos_terminal != 0:
            equations[cccs.pos_terminal] -= i
        if cccs.neg_terminal != 0:
            equations[cccs.neg_terminal] += i

        passing_component

    for vccs in Vccss:
        V_p_node = nodes_Voltages_sym[vccs.pos_terminal]
        V_n_node = nodes_Voltages_sym[vccs.neg_terminal]
        pos_terminal
        neg_terminal
        control_pos
        control_neg
        coefficient
        passing_component
    unknowns = []
    for s in nodes_Voltages_sym:
        unknowns.append(s)
    for s in current_in_vs_sym:
        unknowns.append(s)

    # print("equations:", equations)
    # print("unknowns", unknowns)
    sol = list(solve(equations, unknowns).items())
    # print(cmath.polar(sol[2][1], ))
    k = 0
    while k <= no_nodes:
        value = sol[k][1]
        polar_val = cmath.polar(value)

        print("V at node " + str(k), " = ", round(sol[k][1], 3), " or in Polar (", round(polar_val[0], 3), ",",
              round(polar_val[1] * 180 / cmath.pi, 3), ")", sep="")
        k += 1

    f.close()

    # A, B = linear_eq_to_matrix(equations, unknowns)
    # print(A,B)
    # print((A ** -1)*B)

    # solve_linear_system_LU(equations,unknowns)
