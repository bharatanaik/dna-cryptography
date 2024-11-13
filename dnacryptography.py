import ast ,binascii ,math ,random ,re ,string

class Utils:
    key_del = "<key>"
    no_rounds_del = "<no_rounds>"
    round_del = "<round>"
    reshape_del = "<reshape>"
    crossover_del = "<crossover>"
    crossover_type_del = "<type>"
    single_point_crossover_del = "<single_point>"
    rotate_crossover_del = "<rotate>"
    rotation_offset_del = "<rotation_offset>"
    rotation_types_del = "<rotation_types>"
    mutation_del = "<mutation>"
    complement_mutation_del = "<complement_mutation>"
    alter_mutation_del = "<alter_mutation>"
    mutation_table_del = "<mutation_table>"
    chromosome_del = "<chromosome>"
    # generate encoding tables domains
    two_bit_list = ['00', '01', '10', '11']
    dna_bases = ['A', 'C', 'G', 'T']
    four_bit_list = ['0000', '0001', '0010', '0011', '0100', '0101', '0110', '0111', '1000', '1001', '1010', '1011', '1100', '1101', '1110', '1111']
    two_dna_bases = ['TA', 'TC', 'TG', 'TT', 'GA', 'GC', 'GG', 'GT', 'CA', 'CC', 'CG', 'CT', 'AA', 'AC', 'AG', 'AT']
    # encoding tables and their reversal

    def str2bin(self, sstring):
        return "".join([bin(ord(c))[2:].zfill(8) for c in sstring])

    def bin2str(self, bs):
        return binascii.unhexlify('%x' % int(bs, 2))

    def byte2bin(self, byte_val):
        return bin(byte_val)[2:].zfill(8)

    def bitxor(self, a, b):
        return "".join([str(int(x) ^ int(y)) for (x, y) in zip(a, b)])
    
    def divisors(self, n):
        divs = []
        for i in range(2, int(math.sqrt(n)) + 1):
            if n % i == 0:
                divs.extend([i, n / i])
        divs = [int(d) for d in divs]
        return list(set(divs))
    
    def generate_pre_processing_tables(self):
        self.two_bits_to_dna_base_table = dict(zip(self.two_bit_list, self.dna_bases))
        self.dna_base_to_two_bits_table = dict(zip(self.two_bits_to_dna_base_table.values(), self.two_bits_to_dna_base_table.keys()))

    def generate_mutation_tables(self):
        self.four_bits_to_two_dna_base_table = dict(zip(self.four_bit_list, self.two_dna_bases))
        self.two_dna_base_to_four_bits_table = dict(zip(self.four_bits_to_two_dna_base_table.values(), self.four_bits_to_two_dna_base_table.keys()))

    def group_bits(self, byte, step=2):
        return [byte[i:i + step] for i in range(0, len(byte), step)]
    
    def group_bases(self, dna_seq, step=2):
        bases_groups = []
        for i in range(0, len(dna_seq), step):
            bases_groups.append(dna_seq[i:i + step])
        return bases_groups

    def generate_bits(self, byte_data):
        """
        Take every byte for sequence and group its bits
        :return:
        """
        grouped_bits_data = []
        for byte in byte_data:
            grouped_bits_data.extend(self.group_bits(byte))
        return grouped_bits_data
    
    def binarized_data(self, data):
        # convert every char to ASCII and then to binary
        byte_data = [self.byte2bin(ord(c)) for c in data]
        return self.generate_bits(byte_data)

    def bits_to_dna(self, data, conversion_table):
        # convert binary sequence to DNA sequence
        return "".join([conversion_table[bits] for bits in data])

    def dna_to_bits(self, data, conversion_table):
        # convert DNA sequence to binary sequence
        return "".join([conversion_table[dna_base] for dna_base in data])

    def get_pattern(self, delimiter, s):
        """
        Get the pattern info between delimiters from the string
        """
        regex = "%s(.*?)%s" % (delimiter, delimiter)
        return re.findall(regex, s)

    def encrypt_key(self, data, key):
        """
        Encrypt data with key: data XOR key.
        """
        # repeat key ONLY if data is longer than key and encrypt
        if len(data) > len(key):
            factor = int(len(data) / len(key))
            key += key * factor
            return self.bitxor(data, key)
        return self.bitxor(data, key)
    
    def complement(self, chromosome, point1, point2):
        """
        Flip chromosome bits between point1 and point2.
        """
        new_chromosome = ""

        for i in range(len(chromosome)):
            if i >= point1 and i <= point2:
                if chromosome[i] == '0':
                    new_chromosome += '1'
                else:
                    new_chromosome += '0'
            else:
                new_chromosome += chromosome[i]

        return new_chromosome
    
    def reverse_reshape(self, population):
        # convert the chromosome population back to DNA sequence
        return "".join(population)

class Encryption(Utils):
    def __init__(self)->None:
        self.rounds_no = random.randrange(3, 12, 2)
        self.decryption_key = ""

    def reshape(self, dna_sequence):
        """
        Generate chromosome population.
        :param dna_sequence: a string sequence of DNA bases
        :return: an array of chromosomes, chromosome population
        """
        # choose population size and chromosome length
        divs = self.divisors(len(dna_sequence))
        chromosome_no = divs[random.randint(0, len(divs) - 1)]
        self.chromosome_length = int(len(dna_sequence) / chromosome_no)
        chromosomes = []
        self.decryption_key += self.reshape_del + str(self.chromosome_length) + self.reshape_del
        # retrieve the population
        for i in range(0, len(dna_sequence), self.chromosome_length):
            chromosomes.append(dna_sequence[i:i + self.chromosome_length])
        return chromosomes

    def rotate_crossover(self, population):
        """
        Rotate every chromosome in population left / right according to probability p.
        """
        new_population = []
        self.decryption_key += self.rotate_crossover_del
        # predefined rotation value, varied every round
        rotation_offset = random.randint(1, self.chromosome_length)
        self.decryption_key += self.rotation_offset_del + str(rotation_offset) + self.rotation_offset_del
        self.decryption_key += self.rotation_types_del

        for chromosome in population:
            p = random.uniform(0, 1)
            if p > 0.5:
                self.decryption_key += "right|"
                right_first = chromosome[0: len(chromosome) - rotation_offset]
                right_second = chromosome[len(chromosome) - rotation_offset:]
                new_population.append(right_second + right_first)
            else:
                self.decryption_key += "left|"
                left_first = chromosome[0: rotation_offset]
                left_second = chromosome[rotation_offset:]
                new_population.append(left_second + left_first)

        self.decryption_key += self.rotation_types_del
        self.decryption_key += self.rotate_crossover_del
        return new_population
    
    def single_point_crossover(self, population):
        """
        Combine each two chromosomes in population by using single point crossover.
        """
        self.decryption_key += self.single_point_crossover_del
        new_population = []
        for i in range(0, len(population) - 1, 2):
            candidate1 = population[i]
            candidate2 = population[i + 1]
            length = len(candidate1)
            crossover_point = random.randint(0, length - 1)
            self.decryption_key += str(crossover_point) + "|"
            offspring1 = candidate2[0: crossover_point] + candidate1[crossover_point:]
            offspring2 = candidate1[0: crossover_point] + candidate2[crossover_point:]
            new_population.append(offspring1)
            new_population.append(offspring2)
        # append last chromosome if odd population size
        if len(population) % 2 == 1:
            new_population.append(population[len(population) - 1])
        self.decryption_key += self.single_point_crossover_del
        return new_population
    
    def crossover(self, population):
        p = random.uniform(0, 1)

        if p < 0.33:
            self.decryption_key += self.crossover_type_del + "rotate_crossover" + self.crossover_type_del
            return self.rotate_crossover(population)
        elif p >= 0.33 and p < 0.66:
            self.decryption_key += self.crossover_type_del + "single_point_crossover" + self.crossover_type_del
            return self.single_point_crossover(population)
        else:
            self.decryption_key += self.crossover_type_del + "both" + self.crossover_type_del
            population = self.rotate_crossover(population)
            return self.single_point_crossover(population)
 
    def alter_dna_bases(self, bases:list):
        """
        Alter DNA bases to another one randomly.(e.g. C->G and A->T and viceversa)
        """
        alter_dna_table = {}

        for _ in range(2):
            # choose one randomly then remove it from list
            base1 = bases[random.randint(0, len(bases) - 1)]
            bases.remove(base1)

            # choose one randomly then remove it from list
            base2 = bases[random.randint(0, len(bases) - 1)]
            bases.remove(base2)

            # assign the first to the other
            alter_dna_table[base1] = base2
            alter_dna_table[base2] = base1

        return alter_dna_table

    def mutation(self, population):
        """
        Apply mutation operator by using "complement" and "alter_dna_bases"
        """

        bases = ['A', 'C', 'G', 'T']
        alter_dna_table = self.alter_dna_bases(bases)

        self.decryption_key += self.mutation_table_del + str(alter_dna_table) + self.mutation_table_del

        new_population = []
        for chromosome in population:
            self.decryption_key += self.chromosome_del

            # apply the complement
            b_chromosome = self.dna_to_bits(chromosome, self.dna_base_to_two_bits_table)
            self.decryption_key += self.complement_mutation_del
            point1 = random.randint(0, len(b_chromosome) - 1)
            point2 = random.randint(point1, len(b_chromosome) - 1)
            self.decryption_key += "(%s, %s)" % (point1, point2)
            self.decryption_key += self.complement_mutation_del
            b_chromosome = self.complement(b_chromosome, point1, point2)

            # convert each 4 bits in chromosome to two dna bases using four_bits_to_two_dna_base_table
            four_bits_vector = self.group_bits(b_chromosome, 4)

            last_dna_base = None
            # if the last element is of length 2, don't convert it
            if len(four_bits_vector[len(four_bits_vector) - 1]) == 2:
                last_dna_base = self.two_bits_to_dna_base_table[four_bits_vector[len(four_bits_vector) - 1]]
                # convert only the 4 bits elements
                four_bits_vector = four_bits_vector[:-1]

            dna_seq = self.bits_to_dna(four_bits_vector, self.four_bits_to_two_dna_base_table)
            if last_dna_base is not None:
                dna_seq += last_dna_base
            # and then alter the dna bases between point1 and point2
            self.decryption_key += self.alter_mutation_del
            point1 = random.randint(0, len(dna_seq) - 1)
            point2 = random.randint(point1, len(dna_seq) - 1)
            self.decryption_key += "(%s, %s)" % (point1, point2)
            self.decryption_key += self.alter_mutation_del
            new_chromosome = ""
            for i in range(len(dna_seq)):
                if i >= point1 and i <= point2:
                    new_chromosome += alter_dna_table[dna_seq[i]]
                else:
                    new_chromosome += dna_seq[i]

            new_population.append(new_chromosome)

            self.decryption_key += self.chromosome_del

        return new_population
    
    def dna_get(self, text, key):
        b_data1 = self.binarized_data(text)
        dna_seq = self.bits_to_dna(b_data1, self.two_bits_to_dna_base_table)
        b_data2 = dna_seq
        self.decryption_key += self.no_rounds_del + str(self.rounds_no) + self.no_rounds_del
        # run the algorithm "rounds_no" times
        while self.rounds_no > 0:
            self.decryption_key += self.round_del
            # encrypt data with key after reshaping it back to binary sequence and then convert it back to dna sequence
            b_data2 = self.bits_to_dna(
                self.group_bits(self.encrypt_key(self.dna_to_bits(self.reverse_reshape(b_data2), self.dna_base_to_two_bits_table), key)),
                self.two_bits_to_dna_base_table)
            # create the chromosome population
            b_data2 = self.reshape(b_data2)
            # apply crossover on population
            self.decryption_key += self.crossover_del
            b_data2 = self.crossover(b_data2)
            self.decryption_key += self.crossover_del
            # apply mutation on population
            self.decryption_key += self.mutation_del
            b_data2 = self.mutation(b_data2)
            self.decryption_key += self.mutation_del
            self.rounds_no -= 1
            self.decryption_key += self.round_del
        return self.reverse_reshape(b_data2)
    
    def encrypt(self, text:str):
        key = self.str2bin(''.join(random.SystemRandom().choice(string.ascii_letters + string.digits) for _ in range(16)))

        # append the key first
        self.decryption_key += self.key_del + key + self.key_del
        # generate the encoding tables
        self.generate_pre_processing_tables()
        self.generate_mutation_tables()
        encrypted_text = self.dna_get(text, key)
        return encrypted_text, self.decryption_key


class Decryption(Utils):

    def reshape(self, dna_sequence, reshape_info):
        """
        Generate chromosome population.
        :param dna_sequence: a string sequence of DNA bases
        :param reshape_info: the length of each chromosome
        :return: an array of chromosomes, chromosome population
        """
        chromosome_length = int(reshape_info[0])
        chromosomes = []
        # retrieve the population
        for i in range(0, len(dna_sequence), chromosome_length):
            chromosomes.append(dna_sequence[i:i + chromosome_length])
        return chromosomes

    def rotate_crossover(self, population, rotate_info):
        """
        Rotate every chromosome in population left / right according to probability p.
        """
        new_population = []
        # get the rotation value
        rotation_offset = int(self.get_pattern(self.rotation_offset_del, rotate_info)[0])
        rotations = self.get_pattern(self.rotation_types_del, rotate_info)[0].split("|")[:-1]

        for i in range(len(population)):
            chromosome = population[i]
            direction = rotations[i]

            if direction == "left":
                right_first = chromosome[0: len(chromosome) - rotation_offset]
                right_second = chromosome[len(chromosome) - rotation_offset:]
                new_population.append(right_second + right_first)
            elif direction == "right":
                left_first = chromosome[0: rotation_offset]
                left_second = chromosome[rotation_offset:]
                new_population.append(left_second + left_first)
        return new_population

    def single_point_crossover(self, population, single_point_info):
        """
        Combine each two chromosomes in population by using single point crossover.
        """
        crossover_points = [int(p) for p in single_point_info.split("|") if p != '']
        new_population = []
        for i in range(0, len(population) - 1, 2):
            candidate1 = population[i]
            candidate2 = population[i + 1]
            # get the crossover_point
            crossover_point = crossover_points[int(i / 2)]
            offspring1 = candidate2[0: crossover_point] + candidate1[crossover_point:]
            offspring2 = candidate1[0: crossover_point] + candidate2[crossover_point:]
            new_population.append(offspring1)
            new_population.append(offspring2)

        # append last chromosome if odd population size
        if len(population) % 2 == 1:
            new_population.append(population[len(population) - 1])

        return new_population

    def crossover(self, population, crossover_info):
        # get the crossover type
        crossover_type = self.get_pattern(self.crossover_type_del, crossover_info)[0]

        if crossover_type == "rotate_crossover":
            rotate_info = self.get_pattern(self.rotate_crossover_del, crossover_info)[0]
            return self.rotate_crossover(population, rotate_info)
        elif crossover_type == "single_point_crossover":
            single_point_info = self.get_pattern(self.single_point_crossover_del, crossover_info)[0]
            return self.single_point_crossover(population, single_point_info)
        elif crossover_type == "both":
            rotate_info = self.get_pattern(self.rotate_crossover_del, crossover_info)[0]
            single_point_info = self.get_pattern(self.single_point_crossover_del, crossover_info)[0]
            population = self.single_point_crossover(population, single_point_info)
            return self.rotate_crossover(population, rotate_info)

    def mutation(self, population, mutation_info):
        """
        Apply mutation operator by using "complement" and "alter_dna_bases"
        """

        # extract the alteration table
        alter_dna_table = ast.literal_eval(self.get_pattern(self.mutation_table_del, mutation_info[0])[0])

        chromosomes_info = self.get_pattern(self.chromosome_del, mutation_info[0])

        new_population = []
        for i in range(len(population)):
            chromosome = population[i]
            chromosome_info = chromosomes_info[i]

            # alter back the dna bases between point1 and point2
            alter_info = self.get_pattern(self.alter_mutation_del, chromosome_info)[0]
            point1, point2 = ast.literal_eval(alter_info)
            new_chromosome = ""
            for i in range(len(chromosome)):
                if i >= point1 and i <= point2:
                    new_chromosome += alter_dna_table[chromosome[i]]
                else:
                    new_chromosome += chromosome[i]

            two_bases_vector = self.group_bases(new_chromosome)

            # last base was not converted using four_bits_to_two_dna_base_table
            # convert it to bits using dna_base_to_two_bits_table
            last_two_bits = None
            if len(new_chromosome) % 2 == 1:
                last_two_bits = self.dna_base_to_two_bits_table[new_chromosome[-1]]

                two_bases_vector = two_bases_vector[:-1]

            bits_seq = self.dna_to_bits(two_bases_vector, self.two_dna_base_to_four_bits_table)

            if last_two_bits is not None:
                bits_seq += last_two_bits

            complement_info = self.get_pattern(self.complement_mutation_del, chromosome_info)[0]
            point1, point2 = ast.literal_eval(complement_info)
            b_chromosome = self.complement(bits_seq, point1, point2)
            b_chromosome = self.group_bits(b_chromosome)
            new_chromosome = self.bits_to_dna(b_chromosome, self.two_bits_to_dna_base_table)

            new_population.append(new_chromosome)

        return new_population

    def dna_gdt(self, text, key):

        rounds_no = int(self.get_pattern(self.no_rounds_del, key)[0])
        rounds = self.get_pattern(self.round_del, key)

        b_data = text
        
        # run the algorithm "rounds_no" times
        while rounds_no > 0:
            round_info = rounds[rounds_no - 1]
            # create the chromosome population
            b_data = self.reshape(b_data, self.get_pattern(self.reshape_del, round_info))
            # apply mutation on population
            b_data = self.mutation(b_data, self.get_pattern(self.mutation_del, round_info))
            # apply crossover on population
            b_data = self.crossover(b_data, round_info)
            # decrypt data with key after reshaping it back to binary sequence and then convert it back to dna sequence
            # where decrypt = encrypt(encrypt(data, key), key) and encrypt => xor operation, because (a xor b) xor b = a
            encryption_key = self.get_pattern(self.key_del, key)[0]
            b_data = self.bits_to_dna(
                        self.group_bits(
                            self.encrypt_key(
                                self.dna_to_bits(
                                    self.reverse_reshape(b_data), 
                                    self.dna_base_to_two_bits_table
                                ), 
                                encryption_key
                            )
                        ),
                        self.two_bits_to_dna_base_table
                    )
            rounds_no -= 1

        return self.bin2str(self.dna_to_bits(b_data, self.dna_base_to_two_bits_table)).decode()

    def decrypt(self, encrypted_text:str, key:str):
        self.generate_pre_processing_tables()
        self.generate_mutation_tables()
        decrypted_text = self.dna_gdt(encrypted_text, key)
        return decrypted_text
 