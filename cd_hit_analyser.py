#Importar Módulos
##Para colocar flags em linha de comando
import argparse

##Para executar comandos em bash
import os

##Para utilizar expressoes regulares, para encontrar corrêspondência de identificador em strings
import re

#Criando parser de argumentos (Define, processa e armazena argumentos)
parser = argparse.ArgumentParser()

#Adicionando argumentos ao parser
parser.add_argument("-i", "--input_file", help = "Arquivo Fasta.", required = True)
parser.add_argument("-o", "--output_file", help = "Nome de output .cdhit.", required = True)
parser.add_argument("-id", "--identity", help = "Grau de identidade.", required = True, type = float)
parser.add_argument("-c1s", "--cluster_1_seq", help = "y/n para contar e armazenar clusters com apenas 1 sequência.", choices=["y", "n"], required = False)
parser.add_argument("-t", "--seq_type", help = "Tipo de sequências, 'a' para aminoácidos, 'n' para nucleotídeos. ", choices=["a", "n"], required = True)

#Argumentos para busca de identificadores (Máximo de 5 identificadores)
parser.add_argument("-ident_1", "--identificador_1", help = "Identificador a ser procurado no header das sequências (Pelo menos 2 identificadores são necessários)", required = False)
parser.add_argument("-ident_2", "--identificador_2", help = "Identificador a ser procurado no header das sequências", required = False)
parser.add_argument("-ident_3", "--identificador_3", help = "Identificador a ser procurado no header das sequências", required = False)
parser.add_argument("-ident_4", "--identificador_4", help = "Identificador a ser procurado no header das sequências", required = False)
parser.add_argument("-ident_5", "--identificador_5", help = "Identificador a ser procurado no header das sequências", required = False)


#Analisando argumentos fornecidos na linha de comando e e atribui como objeto args
args = parser.parse_args()

#Salvando argumento em variáveis
arquivo_fasta = args.input_file 
identidade = args.identity
nome_output = args.output_file
cluster_1_seq_y_n = args.cluster_1_seq #Se não houver fornecido, será salvo como none, não dando erro
tipo_seq = args.seq_type

#Se apenas 1 argumento for fornecido, dará erro
if args.identificador_1 and not args.identificador_2:
	parser.error("--identificador_2 é obrigatório caso --identificador_1 seja fornecido.")

#Salvando identificadores em variáveis
id_1 = args.identificador_1
id_2 = args.identificador_2
id_3 = args.identificador_3
id_4 = args.identificador_4
id_5 = args.identificador_5

#Rodar cd-hit
##Fazer correspondência adequada de tamanho de fonte
if tipo_seq == 'n':
	### Se for menor que 1 e maior ou igual a 0.9 tamanho de fonte deve ser 9 
	if 0.9 <= identidade < 1.0:
		os.system(f"cd-hit-est -i {arquivo_fasta} -o {nome_output} -d 0 -T 16 -g 0 -M 75000 -aL 0.97 -aS 0.97 -c {identidade} -n 9 -b 1")

	### Se for menor que 0.9 e maior ou igual a 0.88 tamanho de fonte deve ser 7
	elif 0.88 < identidade < 0.9:
		os.system(f"cd-hit-est -i {arquivo_fasta} -o {nome_output} -d 0 -T 16 -g 0 -M 75000 -aL 0.97 -aS 0.97 -c {identidade} -n 7 -b 1")

	### Se for menor que 0.88 e maior ou igual a 0.85 tamanho de fonte deve ser 6
	elif 0.85 <= identidade < 0.88:
		os.system(f"cd-hit-est -i {arquivo_fasta} -o {nome_output} -d 0 -T 16 -g 0 -M 75000 -aL 0.97 -aS 0.97 -c {identidade} -n 6 -b 1")

	### Se for menor que 0.85 e maior ou igual a 0.8 tamanho de fonte deve ser 5
	elif 0.8 <= identidade < 0.85:
		os.system(f"cd-hit-est -i {arquivo_fasta} -o {nome_output} -d 0 -T 16 -g 0 -M 75000 -aL 0.97 -aS 0.97 -c {identidade} -n 5 -b 1")

        ### Se não tiver em intervalos acima
	else:
		print(f"Erro: Grau de Identidade Inválido \n OBS: Para nucleotídeos o CD-HIT-EST não suporta grau de identidade menor que 0.8")

elif tipo_seq == 'a':
        ### Se for menor que 1 e maior ou igual a 0.7 tamanho de fonte deve ser 5
	if 0.7 <= identidade < 1.0:
		os.system(f"cd-hit -i {arquivo_fasta} -o {nome_output} -d 0 -T 16 -g 0 -M 75000 -aL 0.97 -aS 0.97 -c {identidade} -n 5 -b 1")

	### Se for menor que 0.7 e maior ou igual a 0.6 tamanho de fonte deve ser 4
	elif 0.6 <= identidade < 0.7:
		os.system(f"cd-hit -i {arquivo_fasta} -o {nome_output} -d 0 -T 16 -g 0 -M 75000 -aL 0.97 -aS 0.97 -c {identidade} -n 4 -b 1")

	### Se for menor que 0.6 e maior ou igual a 0.5 tamanho de fonte deve ser 3
	elif 0.5 <= identidade < 0.6:
		os.system(f"cd-hit -i {arquivo_fasta} -o {nome_output} -d 0 -T 16 -g 0 -M 75000 -aL 0.97 -aS 0.97 -c {identidade} -n 3 -b 1")

	### Se for menor que 0.5 e maior ou igual a 0.4 tamanho de fonte deve ser 2
	elif 0.4 <= identidade < 0.5:
		os.system(f"cd-hit -i {arquivo_fasta} -o {nome_output} -d 0 -T 16 -g 0 -M 75000 -aL 0.97 -aS 0.97 -c {identidade} -n 2 -b 1")
        	### Se não tiver em intervalos acima
	else:
        	print(f"Erro: Grau de Identidade Inválido \n OBS: Para aminoácidos o CD-HIT não suporta grau de identidade menor que 0.4") 

else:
	print("Erro: Não informou tipo de sequência")

#Abrir arquivo .clstr e verificar quantidade de cluster
clstr_file = nome_output + ".clstr"

##Função para ler .clstr
def read_clstr(clstr_file):
	#Abrir arquvo . clstr
	with open(clstr_file, 'r') as cl_file:
		#Criar variaveis
		##Quantidade de clusters		
		clusters = 0
		#Maior cluster e quantidade de sequencias nele
		largest_cluster = []
		largest_cluster_quant = 0
		#Menor Cluster e quantidade de sequências nele
		smallest_cluster = []
		smallest_cluster_quant = float("inf") #Começamos com o infinito para que independente do tamanho do cluster o primeiro analizado seja o menor
		#Cluster atual para comparações
		current_cluster = None
		current_cluster_quant = 0
		
		#Lista de quantidade de sequências em cada cluster para fazer média
		quant_seq = 0
		quant_seq_list = []


		#Lista de quantidade de sequências e nomes de clusters de apenas 1 sequência
		if cluster_1_seq_y_n == 'y':
			cluster_1_seq_list = []

		#Criar conjunto para armazenar identificadores de cada cluster
		current_cluster_ids = set() #Não permite elementos duplicados, assim contando apenas aparição de id e não de quantidade

		#Criar dicionário para armazenar quantidade de cluster com os identificadores
		dicionario_cluster_id = {}

		##Loop, no qual a variavel line armazena conteúdo de uma linha do arquivo a cada iteração
		for line in cl_file:
			#Função stars with verifica se linha atual começa com ">Cluster" 
			#A cada linha que comece com ">Cluster" adiciona 1 a contagem
			if line.startswith(">Cluster"):
				clusters += 1
							
			#Verifica se cluster atual é maior ou menor que cluster já analizado
				if current_cluster_quant > largest_cluster_quant:
					largest_cluster_quant = current_cluster_quant
					largest_cluster = [current_cluster]
				
				#Caso seja igual adicionar a lista
				elif current_cluster_quant == largest_cluster_quant:
					largest_cluster.append(current_cluster)

				if 0 < current_cluster_quant < smallest_cluster_quant:
					smallest_cluster_quant = current_cluster_quant
					smallest_cluster = [current_cluster]

				elif current_cluster_quant == smallest_cluster_quant:
					smallest_cluster.append(current_cluster)

				if cluster_1_seq_y_n == 'y' and current_cluster_quant == 1:
					cluster_1_seq_list.append(current_cluster)
					
				#Adicionar contagem de sequencias em cluster em lista
				quant_seq_list.append(current_cluster_quant)
				#Resetar quantidade de sequencias no cluster atual
				current_cluster_quant = 0
				 
				#Adicionar cluster e ids dele em dicionário
				if current_cluster is not None:
					dicionario_cluster_id[current_cluster] = current_cluster_ids
				#Extrair número identificador de cluster		
				current_cluster = int(line.split()[1])

				#Zerar quantidade de ids em conjunto
				current_cluster_ids = set()
	
			#Se linha não indicar novo cluster adiciona a contagem de sequencias dentro do cluster atual
			else:
				current_cluster_quant += 1
				
				if id_1 is not None and re.search(id_1, line):
					current_cluster_ids.add(f"{id_1}")

				elif id_2 is not None and re.search(id_2, line):
					current_cluster_ids.add(f"{id_2}")

				elif id_3 is not None and re.search(id_3, line):
					current_cluster_ids.add(f"{id_3}")

				elif id_4 is not None and re.search(id_4, line):
					current_cluster_ids.add(f"{id_4}")

				elif id_5 is not None and re.search(id_5, line):
					current_cluster_ids.add(f"{id_5}")

		##Checar último cluster por que só efetua comparações no loop se encontrar linha que indica novo cluster
		if current_cluster_quant > largest_cluster_quant:
			largest_cluster_quant = current_cluster_quant
			largest_cluster = current_cluster

		elif current_cluster_quant == largest_cluster_quant:
			largest_cluster.append(current_cluster)

		if 0 < current_cluster_quant < smallest_cluster_quant:
			smallest_cluster_quant = current_cluster_quant
			smallest_cluster = current_cluster

		elif current_cluster_quant == smallest_cluster_quant:
			smallest_cluster.append(current_cluster)

		if cluster_1_seq_y_n == 'y' and current_cluster_quant == 1:
			cluster_1_seq_list.append(current_cluster)

		quant_seq_list.append(current_cluster_quant)

		if id_1 is not None and re.search(id_1, line):
			current_cluster_ids.add(f"{id_1}")

		elif id_2 is not None and re.search(id_2, line):
			current_cluster_ids.add(f"{id_2}")

		elif id_3 is not None and re.search(id_3, line):
			current_cluster_ids.add(f"{id_3}")

		elif id_4 is not None and re.search(id_4, line):
			current_cluster_ids.add(f"{id_4}")

		elif id_5 is not None and re.search(id_5, line):
			current_cluster_ids.add(f"{id_5}")
		
		dicionario_cluster_id[current_cluster] = current_cluster_ids

	#Salvar variáveis da função	
		if cluster_1_seq_y_n != 'y' and id_1 is None:
                        return clusters, largest_cluster, largest_cluster_quant, smallest_cluster, smallest_cluster_quant, quant_seq_list

		elif cluster_1_seq_y_n == 'y' and id_1 is None:
			return clusters, largest_cluster, largest_cluster_quant, smallest_cluster, smallest_cluster_quant, quant_seq_list, cluster_1_seq_list

		elif cluster_1_seq_y_n == 'y' and id_1 is not None:
			return clusters, largest_cluster, largest_cluster_quant, smallest_cluster, smallest_cluster_quant, quant_seq_list, cluster_1_seq_list, dicionario_cluster_id
		
		elif cluster_1_seq_y_n != 'y' and id_1 is not None:
			return clusters, largest_cluster, largest_cluster_quant, smallest_cluster, smallest_cluster_quant, quant_seq_list, dicionario_cluster_id


if cluster_1_seq_y_n != 'y' and id_1 is None:
	quant_clusters, maior_cluster, maior_cluster_quant, menor_cluster, menor_cluster_quant, quant_seq_list  = read_clstr(clstr_file)

elif cluster_1_seq_y_n == 'y' and id_1 is None:
	quant_clusters, maior_cluster, maior_cluster_quant, menor_cluster, menor_cluster_quant, quant_seq_list, cluster_1_seq_list  = read_clstr(clstr_file)

elif cluster_1_seq_y_n == 'y' and id_1 is not None:
	quant_clusters, maior_cluster, maior_cluster_quant, menor_cluster, menor_cluster_quant, quant_seq_list, cluster_1_seq_list, dicionario_cluster_id  = read_clstr(clstr_file)

elif cluster_1_seq_y_n != 'y' and id_1 is not None:
	quant_clusters, maior_cluster, maior_cluster_quant, menor_cluster, menor_cluster_quant, quant_seq_list, dicionario_cluster_id  = read_clstr(clstr_file)

#Calcular média
media = sum(quant_seq_list)/int(quant_clusters)
#Calcular Mediana
##Ordenar lista
quant_seq_list.sort()
##Verificar se quantidade de clusters é impar ou par
if len(quant_seq_list) % 2 !=0: #Impar quando o resto de divisão por 2 é diferente de zero
	mediana = quant_seq_list[len(quant_seq_list) // 2] #// Dá como saída o quociente inteiro da divisão

else: #Par
	elemento1 = quant_seq_list[len(quant_seq_list) // 2 - 1]
	elemento2 = quant_seq_list[len(quant_seq_list) // 2]
	mediana = (elemento1 + elemento2) / 2

#Quantidade de sequencias em menor e menor cluster
largest = len(maior_cluster)
smallest = len(menor_cluster)

#Criando diretório para armazenar saídas
os.makedirs("info_cd_hit", exist_ok= True)

path = "info_cd_hit/"

#Salvar relações calculadas em um arquivo txt
with open (path + "relacoes_cd_hit_" + nome_output + ".txt", "w") as arquivo:
	#Printar grau de identidade
	arquivo.write(f"Grau de identidade fornecido é de {(identidade) * 100}%\n\n")

	arquivo.write(f"Número de clusters encontrados: {quant_clusters}\n\n")
	if len(maior_cluster) > 1 and maior_cluster_quant > 1:
		arquivo.write(f"Maiores clusters | Quantidade de clusters: {largest} | {maior_cluster_quant} sequências\n\n")

	elif len(maior_cluster) > 1 and maior_cluster_quant == 1:
        	arquivo.write(f"Maiores Clusters | Quantidade de clusters: {largest} | {maior_cluster_quant} sequência\n\n")

	elif len(maior_cluster) == 1 and maior_cluster_quant > 1:
		arquivo.write(f"Maior cluster | Quantidade de clusters: {largest} | {maior_cluster_quant} sequências\n\n")

	elif len(maior_cluster) == 1 and maior_cluster_quant == 1:
        	arquivo.write(f"Maior cluster | Quantidade de clusters: {largest} | {maior_cluster_quant} sequência\n\n")

	else:
		print("Erro na análise de quantidade de clusters")

	if len(menor_cluster) > 1 and menor_cluster_quant > 1:
		arquivo.write(f"Menores clusters | Quantidade de clusters: {smallest} | {menor_cluster_quant} sequências\n\n")

	elif len(menor_cluster) > 1 and menor_cluster_quant == 1:
        	arquivo.write(f"Menores clusters | Quantidade de clusters: {smallest} | {menor_cluster_quant} sequência\n\n")

	elif len(menor_cluster) ==1 and menor_cluster_quant > 1:
		arquivo.write(f"Menor cluster | Quantidade de clusters: {smallest} | {menor_cluster_quant} sequências\n\n")

	elif len(menor_cluster) ==1 and menor_cluster_quant == 1:
        	arquivo.write(f"Menor cluster | Quantidade de clusters: {smallest} | {menor_cluster_quant} sequência\n\n")

	else:
        	print("Erro na análise de quantidade de clusters")


	arquivo.write(f"Média da quantidade de sequências por cluster: {media}\n\n")
	arquivo.write(f"Mediana da quantidade de sequências por cluster: {mediana}\n\n")

if cluster_1_seq_y_n == 'y':
	with open (path + "informacao_cluster_1_seq_" + nome_output + ".txt", "w") as arquivo:
        	#Salvar quatidade de clusters com 1 seq e quais são:
		for cluster in cluster_1_seq_list:
			arquivo.write(f">{cluster}\n")


#Salvar Informação de identificadores presentes em um cluster
if id_1 is not None:
	with open (path + "Compartilhamento_clusters_" + nome_output + ".txt", "w") as arquivo:
		for cluster, ids in dicionario_cluster_id.items():
			arquivo.write(f">{cluster}\t{', '.join(ids)}\n") #{', '.join(ids)} para que elementos sejam combinados e separados por linha, removendo elementos adicionais

os.system(f"cat {nome_output}.clstr | awk ' /Cluster/ {{ no=$2;}}; !/Cluster/ {{ id=substr($3, 2, length($3)-4); printf(\"%s\\t%s\\n\", no, id) }} ' > {nome_output}.txt")
