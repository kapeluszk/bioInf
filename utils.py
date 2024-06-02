import heapq
class Graph:

    def __init__(self, max_overlap):
        self.vertices = []
        self.max_overlap = max_overlap
        self.first_vertex = None

    def add_vertex(self, new_vertex):
        self.vertices.append(new_vertex)

    def build(self, spectrum, first_oligonucleotide):
        for oligonucleotide in spectrum:
            new_vertex = Vertex(oligonucleotide)
            self.add_vertex(new_vertex)
            if new_vertex.oligonucleotide == first_oligonucleotide:
                self.first_vertex = new_vertex

        for overlap in range(1, self.max_overlap + 1):
            for source_vertex in self.vertices:
                for target_vertex in self.vertices:
                    if source_vertex.oligonucleotide[overlap:] == target_vertex.oligonucleotide[:len(source_vertex.oligonucleotide) - overlap]:
                        source_vertex.edges.append([target_vertex, overlap, 0])

        print("Graph built correctly...")

class Vertex:
    def __init__(self, oligonucleotide):
        self.oligonucleotide = oligonucleotide
        self.edges = []
        self.visits = 0
        self.attempts = 0

    def increment_visits(self):
        self.visits += 1

    def __lt__(self, other):
        return True

def find_unvisited_vertices(graph, excluded_vertices):
    for vertex in graph.vertices:
        if vertex not in excluded_vertices and vertex.visits == 0:
            if all(adjacent_vertex.visits == 0 for adjacent_vertex, _, _ in vertex.edges):
                if all(adjacent_vertex_2.visits == 0 for adjacent_vertex, _, _ in vertex.edges for adjacent_vertex_2, _, _ in adjacent_vertex.edges):
                    return vertex, False
    return None, len(graph.vertices) == len(excluded_vertices)

def dijkstra_shortest_path(graph, start_vertex, end_vertex):
    distances = {vertex: float('inf') for vertex in graph.vertices}
    distances[start_vertex] = 0
    priority_queue = [(0, start_vertex)]
    parent = {start_vertex: None}

    while priority_queue:
        distance_so_far, current_node = heapq.heappop(priority_queue)
        if current_node == end_vertex:
            break
        for adjacent_vertex, edge_weight, _ in current_node.edges:
            distance = distance_so_far + edge_weight
            if distance < distances[adjacent_vertex]:
                distances[adjacent_vertex] = distance
                parent[adjacent_vertex] = current_node
                heapq.heappush(priority_queue, (distance, adjacent_vertex))

    if distances[end_vertex] == float('inf'):
        return None, float('inf')

    shortest_path = []
    current_node = end_vertex
    while current_node is not None:
        shortest_path.append(current_node)
        current_node = parent[current_node]
    shortest_path.reverse()

    final_path = [[shortest_path.pop(0), 0]]
    for vertex in shortest_path:
        vertex.increment_visits()
        distance = next((edge_weight for adjacent_vertex, edge_weight, _ in final_path[-1][0].edges if adjacent_vertex == vertex), None)
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