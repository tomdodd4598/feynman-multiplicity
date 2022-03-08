# -------------------------------#
# ========== SETTINGS ========== #
# -------------------------------#

# Labels for the incoming and outgoing momenta.
momenta_in = ['p1', 'p2']
momenta_out = ['q1', 'q2']

# Info for the types of vertices in the diagrams.
# Each vertex is defined by a tuple of the number of connections and a label.
vertex_types = [(3, 'g'), (4, 'λ')]

# Maximum number of vertices in the diagrams.
perturbation_order = 2

# Only include diagrams which are normal ordered.
normal_ordering = True

# Only include diagrams which are fully connected.
fully_connected = True

# Only include diagrams which have no vacuum pieces.
no_vacuum_pieces = True

# -------------------------------#
# ========== INTERNAL ========== #
# -------------------------------#


class Vertex:
    def __init__(self, total_connections):
        self.total_connections = total_connections
        self.lines = set()

    def fully_connected(self):
        return len(self.lines) >= self.total_connections

    def add_connection(self, line, is_end_vertex):
        if self.fully_connected():
            return 0
        multiplier = self.total_connections - len(self.lines)
        self.lines.add(line)
        return multiplier if is_end_vertex else 1

    def line_sharing_vertices(self):
        vertices = set()
        for line in self.lines:
            vertices.add(line.other_vertex(self))
        return vertices

    def graph_connected_vertices(self):
        iteration = 0
        vertex_lists = ([self], [])

        def vertex_list(flag):
            return vertex_lists[flag ^ ((iteration & 1) == 1)]

        traversed_vertices = set()
        iter_flag = True
        while iter_flag:
            iter_flag = False
            previous_vertices = vertex_list(0)
            current_vertices = vertex_list(1)
            for vertex in previous_vertices:
                traversed_vertices.add(vertex)
                for line_sharing_vertex in vertex.line_sharing_vertices():
                    if line_sharing_vertex not in traversed_vertices:
                        iter_flag = True
                        current_vertices.append(line_sharing_vertex)
            previous_vertices.clear()
            iteration += 1
        return traversed_vertices


class InternalVertex(Vertex):
    def __init__(self, total_connections, label):
        super().__init__(total_connections)
        self.label = label
        self.original_label = label

    def __str__(self):
        return self.label


def make_internal_vertices_distinguishable(vertices):
    label_dict: dict[str, int] = dict()
    for vertex in vertices:
        count = label_dict.get(vertex.label)
        label_dict[vertex.label] = 1 if count is None else 1 + count
    import copy
    copy_dict = copy.deepcopy(label_dict)
    for vertex in vertices:
        count = copy_dict[vertex.label]
        if count > 1:
            old_label = vertex.label
            vertex.label = vertex.label + str(1 + count - label_dict[old_label])
            label_dict[old_label] -= 1


class ExternalVertex(Vertex):
    def __init__(self, total_connections, incoming, momentum):
        super().__init__(total_connections)
        self.incoming = incoming
        self.momentum = momentum

    def __str__(self):
        return ('' if self.incoming else '-') + self.momentum


class InVertex(ExternalVertex):
    def __init__(self, total_connections, momentum):
        super().__init__(total_connections, True, momentum)


class OutVertex(ExternalVertex):
    def __init__(self, total_connections, momentum):
        super().__init__(total_connections, False, momentum)


external_vertices = [InVertex(1, p) for p in momenta_in]
external_vertices += [OutVertex(1, q) for q in momenta_out]


class Line:
    def __init__(self, vertex_start, vertex_end):
        self.vertex_start = vertex_start
        self.vertex_end = vertex_end

    def other_vertex(self, vertex):
        return self.vertex_start if vertex is self.vertex_end else self.vertex_end

    def __str__(self):
        return '{' + str(self.vertex_start) + ',' + str(self.vertex_end) + '}'


class Diagram:
    def __init__(self, internal_vertices):
        make_internal_vertices_distinguishable(internal_vertices)
        self.vertices: list[Vertex] = external_vertices + internal_vertices
        self.vertex_count = len(self.vertices)
        connections = self.total_vertex_connections()
        self.is_valid = (connections & 1) == 0
        self.iters_remaining = connections // 2
        self.lines = []
        self.multiplier = 1

    def total_vertex_connections(self):
        count = 0
        for vertex in self.vertices:
            count += vertex.total_connections
        return count

    def add_line(self, start_index, end_index):
        start_vertex = self.vertices[start_index]
        end_vertex = self.vertices[end_index]
        line = Line(start_vertex, end_vertex)
        self.lines.append(line)
        self.multiplier *= start_vertex.add_connection(line, False)
        self.multiplier *= end_vertex.add_connection(line, True)

    def iterate(self):
        diagrams = []
        unfinished = self.iters_remaining > 0
        if unfinished:
            start_index = -1
            for index in range(self.vertex_count):
                if not self.vertices[index].fully_connected():
                    start_index = index
                    break
            if start_index >= 0:
                for end_index in range(start_index, self.vertex_count):
                    import copy
                    diagram = copy.deepcopy(self)
                    diagram.add_line(start_index, end_index)
                    if diagram.multiplier > 0:
                        diagram.iters_remaining -= 1
                        diagrams.append(diagram)
        else:
            diagrams.append(self)
        return unfinished, diagrams

    def has_one_vertex_loop(self):
        for line in self.lines:
            if line.vertex_start is line.vertex_end:
                return True
        return False

    def has_vacuum_piece(self):
        for vertex in self.vertices:
            external_connection = False
            for connected_vertex in vertex.graph_connected_vertices():
                if isinstance(connected_vertex, ExternalVertex):
                    external_connection = True
                    break
            if not external_connection:
                return True
        return False

    def is_fully_connected(self):
        vertex_count = len(self.vertices)
        if vertex_count == 0:
            return False
        return len(self.vertices[0].graph_connected_vertices()) == vertex_count

    def __str__(self):
        s = str(self.multiplier) if self.multiplier > 1 else ''
        s += '⟨'
        for line in self.lines:
            s += str(line)
        return s + '⟩'


def internal_vertex_perms(order):
    import itertools
    return itertools.combinations_with_replacement(vertex_types, order)


def main():
    iteration = 0
    diagram_lists = ([], [])

    def diagram_list(flag):
        return diagram_lists[flag ^ ((iteration & 1) == 1)]

    initial_diagrams = diagram_list(0)
    for order in range(1 + perturbation_order):
        for vertex_perm in internal_vertex_perms(order):
            diagram = Diagram([InternalVertex(connections, label) for connections, label in vertex_perm])
            if diagram.is_valid:
                initial_diagrams.append(diagram)

    iter_flag = True
    while iter_flag:
        iter_flag = False
        previous_diagrams = diagram_list(0)
        current_diagrams = diagram_list(1)
        for diagram in previous_diagrams:
            unfinished, next_diagrams = diagram.iterate()
            iter_flag |= unfinished
            current_diagrams += next_diagrams
        previous_diagrams.clear()
        iteration += 1

    all_diagrams = diagram_list(0)
    filtered_diagrams = []

    for diagram in all_diagrams:
        valid = True
        if valid and normal_ordering:
            valid = not diagram.has_one_vertex_loop()

        if valid and fully_connected:
            valid = diagram.is_fully_connected()

        if valid and no_vacuum_pieces:
            valid = not diagram.has_vacuum_piece()

        if valid:
            filtered_diagrams.append(diagram)

    for diagram in filtered_diagrams:
        print(diagram)


if __name__ == '__main__':
    main()
