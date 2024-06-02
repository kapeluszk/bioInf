import random
import math
import os


class Instance:
    def __init__(self):
        self.n = None
        self.k = None
        self.negativeAmount = None
        self.positiveAmount = None
        self.spectrum = list()
        self.repeats = 0
        self.dna = ""
        self.first = ""

    def load_instance(self, file_path):
        with open(file_path, "r") as file:
            self.dna = file.readline().rstrip()
            self.first = file.readline().rstrip()
            data = file.readline().split()
            print(data)
            self.n = int(data[0])
            self.k = int(data[1])
            self.repeats = int(data[2])
            self.negative = float(data[3])
            self.negativeAmount = math.floor((self.n - self.k + 1) * self.negative)
            self.positive = float(data[4])
            self.positiveAmount = math.floor((self.n - self.k + 1) * self.positive)
            self.spectrum = [line.rstrip() for line in file]

    def load_DNA_from_txt(self, file_path):
        with open(file_path, "r") as file:
            dna = file.readline().rstrip()
            n = int(input("Enter which first n-elements of DNA you want to use: "))
            self.n = n
            self.dna = dna[:self.n]
            k = input("Enter the length of oligonucleotide: ")
            self.k = int(k)
            self.first = self.dna[:self.k]
            self.spectrum = [self.dna[i:i + self.k] for i in range(n - self.k + 1)]
            self.spectrum = [self.first] + sorted(self.spectrum[1:])
            self.repeats = sum(self.spectrum[i] == self.spectrum[i+1] for i in range(len(self.spectrum)-1))
            self.negative = 0
            self.negativeAmount = math.floor((n - self.k + 1) * self.negative) + self.repeats
            self.negative = self.negativeAmount / (n - self.k + 1)
            self.positiveAmount = 0
            self.positive = 0
            print("DNA: ", self.dna)
            print("First oligonucleotide: ", self.first)
            print("Spectrum: ", self.spectrum)



    def create_instance(self, n, k):
        self.n = n
        self.k = k
        self.positive = 0
        self.positiveAmount = 0
        self.negative = 0
        while True:
            self.dna = ''.join(random.choices('ACGT', k=self.n))
            self.first = self.dna[:self.k]
            self.spectrum = [self.dna[i:i + self.k] for i in range(self.n - self.k + 1)]
            self.spectrum = [self.first] + sorted(self.spectrum[1:])
            self.repeats = sum(self.spectrum[i] == self.spectrum[i + 1] for i in range(len(self.spectrum) - 1))
            if self.repeats == 0:
                break
        print("Repeats: ", self.repeats)
        self.negativeAmount = math.floor((self.n - self.k + 1) * self.negative) + self.repeats
        self.negative = self.negativeAmount / (self.n - self.k + 1)
        print("DNA: ", self.dna)

    def save_instance(self, folder_path, filename):
        if not os.path.exists(folder_path):
            os.makedirs(folder_path)
        file_path = os.path.join(folder_path, filename)
        with open(file_path, "w+") as file:
            file.write(self.dna + "\n")
            file.write(self.first + "\n")
            file.write(f"{self.n} {self.k} {self.repeats} {self.negative} {self.positive}\n")
            for s in self.spectrum:
                file.write(s + "\n")

    def insert_positive_errors(self, amount):
        for _ in range(amount):
            index, position = self._get_random_index_and_position()
            char = self._get_random_char(index, position)
            new_oligonucleotide = self.spectrum[index][:position] + char + self.spectrum[index][position + 1:]
            while new_oligonucleotide in self.spectrum:
                index, position = self._get_random_index_and_position()
                char = self._get_random_char(index, position)
                new_oligonucleotide = self.spectrum[index][:position] + char + self.spectrum[index][position + 1:]
            print(f"Replacing {self.spectrum[index][position]} with {char} at position {position} in {self.spectrum[index]}")
            self.spectrum[index] = new_oligonucleotide
        self.positiveAmount += amount
        self.positive = self.positiveAmount / (self.n - self.k + 1)

    def _get_random_index_and_position(self):
        index = random.randint(3, len(self.spectrum) - 1)
        position = random.randint(0, self.k - 1)
        return index, position

    def _get_random_char(self, index, position):
        return random.choice("ACGT".replace(self.spectrum[index][position], ""))

    def _is_oligonucleotide_in_spectrum(self, index, char, position):
        return self.spectrum[index][:position] + char + self.spectrum[index][position + 1:] in self.spectrum

    def insert_negative_errors(self, amount):
        for _ in range(amount - self.repeats):
            index = random.randint(1, len(self.spectrum) - 1)
            print(f"Deleting {self.spectrum[index]}")
            self.spectrum.pop(index)
        self.negativeAmount += amount
        self.negative = self.negativeAmount / (self.n - self.k + 1)