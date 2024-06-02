import generator
import utils
import os


def reconstruct_DNA(graph, dna_length, instance):
    visited_nodes = {graph.first}
    path = [[graph.first, 0]]
    path_length = len(graph.vertices[0].oligonucleotide)
    end_mode = False
    double_visited_count = 0
    MAX_ATTEMPTS = 3

    while path_length < dna_length:
        if len(visited_nodes) / len(graph.vertices) < 0.8 and not end_mode:
            if path[-1][0].edges:
                checked, path_length, double_visited_count = verify_edges(path, visited_nodes, path_length, double_visited_count, instance, graph, dna_length)
            if not checked:
                path, path_length, end_mode = handle_unchecked(path, visited_nodes, path_length, graph, end_mode)
        else:
            path, path_length = handle_else(path, path_length, end_mode, graph, dna_length)

        print("added vertex: ", path[-1][0].oligonucleotide)
        if path[-1][0].attempts < MAX_ATTEMPTS:
            path[-1][0].attempts += 1
        else:
            #do not pop the first vertex
            if path[-1][0] is not graph.first:
                path.pop()
                path_length -= 1
            if path and path[-1][0] is not None:
                for v in path[-1][0].edges:
                    if v[0] not in [vertex[0] for vertex in path]:
                        if path_length + v[1] <= dna_length:
                            path.append([v[0], v[1]])
                            path_length += v[1]
                            v[0].visits += 1j
                        break

    print("Path length: ", path_length)
    print("Double visited count: ", double_visited_count)
    return path

def verify_edges(path, visited_nodes, path_length, double_visited_count, instance, graph, dna_length):
    checked = False
    if len(path[-1][0].edges) == 0:
        print("No edges")
        path, path_length = handle_no_edges(path, path_length, False, graph, dna_length)
        return checked, path_length, double_visited_count
    for v in path[-1][0].edges:
        if double_visited_count > instance.negativeAmount * 1.5:
            if v[0].visits == 0 and v[1] == 1:
                path, visited_nodes, path_length, double_visited_count = update_path_and_visited_nodes(path, visited_nodes, path_length, v, double_visited_count)
                checked = True
                break
        else:
            if v[1] == 1:
                path, visited_nodes, path_length, double_visited_count = update_path_and_visited_nodes(path, visited_nodes, path_length, v, double_visited_count)
                checked = True
                break
    return checked, path_length, double_visited_count

def update_path_and_visited_nodes(path, visited_nodes, path_length, v, double_visited_count):
    path.append([v[0], 1])
    visited_nodes.add(v[0])
    path_length += v[1]
    v[0].visits += 1
    if v[0].visits > 1:
        double_visited_count += 1
    return path, visited_nodes, path_length, double_visited_count

def handle_unchecked(path, visited_nodes, path_length, graph, end_mode):
    exceptions = []
    exhausted = False
    new_path = None
    while not exhausted:
        dest, exhausted = utils.find_free_spots(graph, exceptions)
        if dest is not None:
            new_path, distance = utils.dijkstra_shortest_path(graph, path[-1][0], dest)
        if new_path is None or dest is None:
            exceptions.append(dest)
        if new_path is not None:
            break

    if new_path is not None:
        path_length += distance
        path = path + new_path[1:]
    if exhausted:
        end_mode = True
    return path, path_length, end_mode

def handle_else(path, path_length, end_mode, graph, dna_length):
    if len(path[-1][0].edges) != 0:
        v = path[-1][0].edges[0]
        path.append([v[0], v[1]])
        path_length += v[1]
        v[0].visits += 1
    else:
        path, path_length = handle_no_edges(path, path_length, end_mode, graph, dna_length)
    return path, path_length

def handle_no_edges(path, path_length, end_mode, graph, dna_length):
    v = path[-1][0]
    found = False
    for i in [4, 5]:
        if found:
            break
        for on in graph.vertices:
            if on.oligonucleotide[i:] == v.oligonucleotide[:len(on.oligonucleotide) - i]:
                path.append([on, i])
                on.visits += 1
                path_length += i
                found = True
                break
    if end_mode and not found:
        path_length = dna_length
    return path, path_length

def path_to_string(path):
    result = ""
    for vertex, weight in path:
        result += ''.join(vertex.oligonucleotide[0])

    result = result[:-1] + path[-1][0].oligonucleotide
    return result

def save_batch_results_to_csv(folder_path, avg_levenstein, levenstein_table):
    # save to a csv in first line avg_levenstein, in second line instance name and levenstein distance
    with open(os.path.join(folder_path, 'results.csv'), 'w') as file:
        file.write(str(avg_levenstein) + '\n')
        for i, levenstein in enumerate(levenstein_table):
            file.write(f'instance{i}.txt, {levenstein}\n')

def main():
    # Create a folder for results
    folder_path = os.path.join("results", "test")
    if not os.path.exists(folder_path):
        os.makedirs(folder_path)

    print('What do you want to do?')
    print('1. Generate an instance and find a path')
    print('2. Load an instance from a file and find a path')
    print('3. Generate an instance, add errors and find a path')
    print('4. Load an instance from a file, add errors and find a path')
    print('5. Generate multiple instances and save them to files')
    print('6. Load multiple instances from files and find paths')
    print('7. Load multiple instances from files, add errors and find paths')
    print('8. Load single line string with real DNA and prepare instance')

    choice = input('Enter the number of the option: ')
    if choice == '1':
        n = int(input('Enter the length of DNA: '))
        k = int(input('Enter the length of oligonucleotide: '))
        # Generate an instance
        instance = generator.Instance()
        instance.create_instance(n, k)
        # Save the instance to a file
        generator.Instance.save_instance(instance, folder_path, 'instance.txt')
        graph = utils.Graph(k-1)
        graph.build(instance.spectrum, instance.first)
        path = reconstruct_DNA(graph, n, instance)
        string_path = path_to_string(path)
        print(string_path)
        print("stringPath length: ", len(string_path))
        dist = utils.levenshteinDistance(string_path, instance.dna)
        print("Levenshtein distance: ", dist)
    elif choice == '2':
        file_path = input('Enter the path to the file: ')
        instance = generator.Instance.load_instance(file_path)
        n = instance.n
        k = instance.k
        print(n,k)
        graph = utils.Graph(k-1)
        graph.build(instance.spectrum, instance.first)
        path = reconstruct_DNA(graph, n, instance)
        string_path = path_to_string(path)
        print(string_path)
        print("stringPath length: ", len(string_path))
        dist = utils.levenshteinDistance(string_path, instance.dna)
        print("Levenshtein distance: ", dist)
    elif choice == '3':
        n = int(input('Enter the length of DNA: '))
        k = int(input('Enter the length of oligonucleotide: '))
        # Generate an instance
        instance = generator.Instance()
        instance.create_instance(n, k)
        # Save the instance to a file
        generator.Instance.save_instance(instance, folder_path, 'instance.txt')
        neg_err = int(input('Enter the number of negative errors: '))
        pos_err = int(input('Enter the number of positive errors: '))

        instance.insert_negative_errors(neg_err)
        instance.insert_positive_errors(pos_err)
        graph = utils.Graph(k-1)
        graph.build(instance.spectrum, instance.first)
        path = reconstruct_DNA(graph, n, instance)
        string_path = path_to_string(path)
        print(string_path)
        print("stringPath length: ", len(string_path))
        dist = utils.levenshteinDistance(string_path, instance.dna)
        print("Levenshtein distance: ", dist)
    elif choice == '4':
        file_path = input('Enter the path to the file: ')
        instance = generator.Instance.load_instance(file_path)
        n = instance.n
        k = instance.k
        neg_err = int(input('Enter the number of negative errors: '))
        pos_err = int(input('Enter the number of positive errors: '))

        instance.insert_negative_errors(neg_err)
        instance.insert_positive_errors(pos_err)
        graph = utils.Graph(k-1)
        graph.build(instance.spectrum, instance.first)
        path = reconstruct_DNA(graph, n, instance)
        string_path = path_to_string(path)
        print(string_path)
        print("stringPath length: ", len(string_path))
        dist = utils.levenshteinDistance(string_path, instance.dna)
        print("Levenshtein distance: ", dist)
    elif choice == '5':
        amount = int(input('Enter the number of instances to generate: '))
        n = int(input('Enter the length of DNA: '))
        k = int(input('Enter the length of oligonucleotide: '))
        for i in range(amount):
            instance = generator.Instance()
            instance.create_instance(n, k)
            generator.Instance.save_instance(instance, folder_path, f'instance{i}.txt')
    elif choice == '6':
        print('files have to have same name with increasing number at the end, starting from 0')
        amount = int(input('Enter the number of instances to load: '))
        file_path = input('Enter the path to the folder: ')
        overlap = int(input('Enter the overlap: '))
        print(file_path)
        avg_levenstein = 0
        lev_table = []
        for i in range(amount):
            instance = generator.Instance()
            instance.load_instance(os.path.join(file_path, f'instance{i}.txt'))
            n = instance.n
            k = instance.k
            print(n, k)
            graph = utils.Graph(overlap)
            graph.build(instance.spectrum, instance.first)
            path = reconstruct_DNA(graph, n, instance)
            string_path = path_to_string(path)
            leven = utils.levenshteinDistance(string_path, instance.dna)
            lev_table.append(leven)
            print("instance", i, "Levenshtein distance: ", leven)
        avg_levenstein = sum(lev_table) / len(lev_table)
        print("Average Levenshtein distance: ", avg_levenstein)
        save_batch_results_to_csv(folder_path, avg_levenstein, lev_table)
    elif choice == '7':
        print('files have to have same name with increasing number at the end, starting from 0')
        amount = int(input('Enter the number of instances to load: '))
        file_path = input('Enter the path to the folder: ')
        neg_err = int(input('Enter the number of negative errors: '))
        pos_err = int(input('Enter the number of positive errors: '))
        avg_levenstein = 0
        lev_table = []
        for i in range(amount):
            instance = generator.Instance()
            instance.load_instance(os.path.join(file_path, f'instance{i}.txt'))
            n = instance.n
            k = instance.k
            print(n, k)
            instance.insert_negative_errors(neg_err)
            instance.insert_positive_errors(pos_err)
            graph = utils.Graph(k-1)
            graph.build(instance.spectrum, instance.first)
            path = reconstruct_DNA(graph, n, instance)
            string_path = path_to_string(path)
            leven = utils.levenshteinDistance(string_path, instance.dna)
            lev_table.append(leven)
            print("instance", i, "Levenshtein distance: ", leven)
        avg_levenstein = sum(lev_table) / len(lev_table)
        print("Average Levenshtein distance: ", avg_levenstein)
    elif choice == '8':
        instance = generator.Instance()
        file_path = input('Enter the path to the file: ')
        instance.load_DNA_from_txt(file_path)
        graph = utils.Graph(instance.k-1)
        graph.build(instance.spectrum, instance.first)
        path = reconstruct_DNA(graph, instance.n, instance)
        string_path = path_to_string(path)
        lev = utils.levenshteinDistance(string_path, instance.dna)
        print("Levenshtein distance: ", lev)
    else:
        print('Invalid choice')

    print("End of the program...")



# Press the green button in the gutter to run the script.
if __name__ == '__main__':
    main()

# See PyCharm help at https://www.jetbrains.com/help/pycharm/
