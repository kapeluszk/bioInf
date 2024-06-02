import heapq
class Graph:

    def __init__(self, max_overlay):
        self.vertices = []
        self.overlay = max_overlay
        self.first = None

    def add_vertex(self, vertex):
        self.vertices.append(vertex)

    def build(self, spectrum, first):
        for oligonucleotide in spectrum:
            v = Vertex(oligonucleotide)
            self.add_vertex(v)
            if v.oligonucleotide == first:
                self.first = v

        for i in range(1, self.overlay + 1):
            for vertex1 in self.vertices:
                for vertex2 in self.vertices:
                    if vertex1.oligonucleotide[i:] == vertex2.oligonucleotide[:len(vertex1.oligonucleotide) - i]:
                        vertex1.edges.append([vertex2, i, 0])

        print("Graph built correctly...")

class Vertex:
    def __init__(self, oligonucleotide):
        self.oligonucleotide = oligonucleotide
        self.edges = []
        self.visits = 0
        self.attempts = 0

    def add_visit(self):
        self.visits += 1

    def __lt__(self, other):
        return True

def find_free_spots(graph, exceptions):
    for vertex in graph.vertices:
        if vertex not in exceptions and vertex.visits == 0:
            if all(adjacent_vertex.visits == 0 for adjacent_vertex, _, _ in vertex.edges):
                if all(adjacent_vertex_2.visits == 0 for adjacent_vertex, _, _ in vertex.edges for adjacent_vertex_2, _, _ in adjacent_vertex.edges):
                    return vertex, False
    return None, len(graph.vertices) == len(exceptions)

def dijkstra_shortest_path(graph, start_vertex, end_vertex):
    distances = {vertex: float('inf') for vertex in graph.vertices}
    distances[start_vertex] = 0
    priority_queue = [(0, start_vertex)]
    parent = {start_vertex: None}

    while priority_queue:
        current_distance, current_vertex = heapq.heappop(priority_queue)
        if current_vertex == end_vertex:
            break
        for neighbor, edge_weight, _ in current_vertex.edges:
            distance = current_distance + edge_weight
            if distance < distances[neighbor]:
                distances[neighbor] = distance
                parent[neighbor] = current_vertex
                heapq.heappush(priority_queue, (distance, neighbor))

    if distances[end_vertex] == float('inf'):
        return None, float('inf')

    shortest_path = []
    current_vertex = end_vertex
    while current_vertex is not None:
        shortest_path.append(current_vertex)
        current_vertex = parent[current_vertex]
    shortest_path.reverse()

    final_path = [[shortest_path.pop(0), 0]]
    for vertex in shortest_path:
        vertex.visits += 1
        distance = next((edge_weight for neighbor, edge_weight, _ in final_path[-1][0].edges if neighbor == vertex), None)
        final_path.append([vertex, distance])
    return final_path, distances[end_vertex]

def levenshteinDistance(str1, str2):
    m = len(str1)
    n = len(str2)
    d = [[i] for i in range(1, m + 1)]   # d matrix rows
    d.insert(0, list(range(0, n + 1)))   # d matrix columns
    for j in range(1, n + 1):
        for i in range(1, m + 1):
            if str1[i - 1] == str2[j - 1]:   # Python (string) is 0-based
                substitutionCost = 0
            else:
                substitutionCost = 1
            d[i].insert(j, min(d[i - 1][j] + 1,
                               d[i][j - 1] + 1,
                               d[i - 1][j - 1] + substitutionCost))
    return d[-1][-1]