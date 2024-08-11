import streamlit as st
from transformers import AutoTokenizer, AutoModelForSequenceClassification
import torch

tokenizer = AutoTokenizer.from_pretrained("ankur0402/dnabert2first",trust_remote_code=True)
model = AutoModelForSequenceClassification.from_pretrained("ankur0402/dnabert2first",trust_remote_code=True)

class MySeq:
    def __init__(self, sequence, alphabet):
        self.sequence = sequence
        self.alphabet = alphabet

    def __str__(self):
        return self.sequence

    def __len__(self):
        return len(self.sequence)

    def get_item(self, index):
        if 0 <= index < len(self.sequence):
            return self.sequence[index]
        else:
            raise IndexError("Index out of range")

    def get_seq_biotype(self):
        return self.alphabet

    def getalphabet(self):
        if self.alphabet == "dna":
            return "ACGT"
        elif self.alphabet == "rna":
            return "ACGU"
        elif self.alphabet == "protein":
            return "ACDEFGHIKLMNPQRSTVWY"
        else:
            return None

class DeterministicMotifFinding:
    """Class for deterministic motif finding."""
    count = 0
    l=[]
    def __init__(self, size=8, seqs=None):
        self.motif_size = size
        self.seqs = seqs if seqs is not None else []
        self.alphabet = None # Initialize alphabet to None
        self.determine_alphabet()

    def __len__(self):
        return len(self.seqs)

    def __getitem__(self, n):
        return self.seqs[n]

    def seq_size(self, i):
        return len(self.seqs[i])

    def read_file(self, fic, t):
        for s in open(fic, "r"):
            self.seqs.append(s.strip().upper())
            print(self.seqs)
        self.determine_alphabet()  # Determine the alphabet after reading sequences from the file

    def determine_alphabet(self):
        # Create a set to store unique characters
        unique_chars = set()
        for seq in self.seqs:
            unique_chars.update(seq)

        # Convert the set of unique characters to a sorted list to define the alphabet
        self.alphabet = ''.join(sorted(unique_chars))

    def create_motif_from_indexes(self, indexes):
        pseqs = []
        res = [[0] * self.motif_size for _ in range(len(self.alphabet))]


        for i, ind in enumerate(indexes):
            if ind + self.motif_size > len(self.seqs[i]):
                continue  # Skip this iteration or handle it as needed
            subseq = self.seqs[i][ind:(ind + self.motif_size)]

            for j in range(self.motif_size):
                for k, char in enumerate(self.alphabet):
                    if subseq[j] == char:
                        res[k][j] += 1


        return res

    def score(self, s):
        """Calculate the score of a motif represented by indices s."""
        score = 0
        mat = self.create_motif_from_indexes(s)
        if not mat or not mat[0]:
            print("Error: Matrix dimensions are not as expected.")
            return  # or handle the error appropriately
        else:
            for j in range(len(mat[0])):
                maxcol = mat[0][j]
                for i in range(1, len(mat)):
                    if mat[i][j] > maxcol:
                        maxcol = mat[i][j]
            score += maxcol

        return score

    def score_multiplicative(self, s):
        """Calculate the multiplicative score of a motif represented by indices s."""
        score = 1.0
        mat = self.create_motif_from_indexes(s)
        for j in range(len(mat[0])):
            maxcol = mat[0][j]
            for i in range(1, len(mat)):
                if mat[i][j] > maxcol:
                    maxcol = mat[i][j]
            score *= maxcol

        return score

    def calculate_alphabet(self):
        unique_chars = set()
        for sequence in self.seqs:
            unique_chars.update(set(sequence))
        return ''.join(sorted(unique_chars))


    def next_solution(self, s):
        next_sol = [0] * len(s)
        pos = len(s) - 1
        while pos >= 0 and s[pos] == self.seq_size(pos) - self.motif_size:
            pos -= 1
        if pos < 0:
            next_sol = None
        else:
            for i in range(pos):
                next_sol[i] = s[i]
            next_sol[pos] = s[pos] + 1
            for i in range(pos + 1, len(s)):
                next_sol[i] = 0
        return next_sol

    def exhaustive_search(self):
        best_score = -1
        res = []
        s = [0] * len(self.seqs)
        while s is not None:
            sc = self.score(s)
            if sc is not None:  # Check if the score is not None
                if sc > best_score:
                    best_score = sc
                    res = s
            s = self.next_solution(s)
        return res


def create_matrix_zeros(nrows, ncols):
    res = []
    for i in range(0, nrows):
        res.append([0] * ncols)
    return res

def print_matrix(mat):
    for i in range(0, len(mat)):
        print(mat[i])

class MyMotifs:
#Class to handle Probabilistic Weighted Matrix
  def __init__( self , seqs = [], pwm = [], alphabet = None):
    if seqs:
      self.size = len (seqs[0])
      self.seqs = seqs # objet from class MySeq
      self.pwm=pwm
      self .alphabet = seqs[0].getalphabet()
      self.do_counts()
      self.create_pwm()
    else :
      self.pwm = pwm
      self.size = len (pwm[0])
      self.alphabet = alphabet


  def __len__ ( self ):
      return self .size

  def do_counts( self ):
      self .counts = create_matrix_zeros( len ( self .alphabet), self .size)
      for s in self.seqs:
        for i in range( self .size):
          char = s.get_item(i)
          lin = self .alphabet.index(char)
          self .counts[lin][i] += 1
  def create_pwm( self ):
      if self .counts == None: self .do_counts()
      self .pwm = create_matrix_zeros( len ( self .alphabet), self .size)
      for i in range( len ( self .alphabet)):
        for j in range( self .size):
          self .pwm[i][j] = float ( self .counts[i][j]) / len ( self .seqs)

  def consensus( self ):
    """ returns the sequence motif obtained with the most
    frequent symbol at each position of the motif"""
    res = ""
    for j in range( self .size):
      maxcol = self .counts[0][j]
      maxcoli = 0
      for i in range (1, len ( self .alphabet) ):
        if self .counts[i][j] > maxcol:
          maxcol = self .counts[i][j]
          maxcoli = i
      res += self .alphabet[maxcoli]
    return res

  def masked_consensus( self ):
    """ returns the sequence motif obtained with the symbol that
    occurs in at least 50% of the input sequences"""
    res = ""
    for j in range( self .size):
      maxcol = self .counts[0][j]
      maxcoli = 0
      for i in range (1, len ( self .alphabet) ):
        if self .counts[i][j] > maxcol:
          maxcol = self .counts[i][j]
          maxcoli = i
      if maxcol > len ( self .seqs) / 2:
        res += self .alphabet[maxcoli]
      else :
        res += "−"
    return res

  def probability_sequence( self , seq):
    res = 1.0
    for i in range( self .size):
      lin = self .alphabet.index(seq[i])
      res *= self .pwm[lin][i]
    return res

  def probability_all_positions( self , seq):  #update to give all probabilities by dividing seqs into multiple seqs of size
    res = []
    for k in range( len (seq)-self .size+1):
      res.append( self .probability_sequence(seq))
    return res
  def most_probable_sequence( self , seq):
    """ Returns the index of the most probable sub−sequence of
    the input sequence"""
    maximum = -1.0
    maxind = -1
    for k in range( len (seq)-self .size):
      p = self .probability_sequence(seq[k:k+ self .size])
      if (p > maximum):
        maximum = p
        maxind = k
    return maxind
  def create_motif( self , seqs):
    l = []
    ind = self .most_probable_sequence(seqs)
    subseq =( seqs[ind:(ind+ self .size)]) #,self.seqs.get_seq_biotype() )
    l.append(subseq)
    return l[0]

def main():
    page = st.sidebar.radio("Select Page", ["Motif Summary", "Motif Finding", "Promoter Prediction"])
    if page == "Motif Summary":
        home_page()
    elif page == "Motif Finding":
        page1()
    elif page == "Promoter Prediction":
        page2()
    # Example usage of your classes
def home_page():
    st.header("Welcome to the Home Page")
    st.write("This is the main page of the app. Use the sidebar to navigate.")
    st.title("Motif Finding Streamlit App")

    # Input for multiple DNA sequences
    st.header("Enter DNA Sequences:")
    sequences = []
    for i in range(1, 9):
        sequence = st.text_input(f"Sequence {i}:", key=f"seq_input_{i}")
        sequences.append(MySeq(sequence, "dna"))

    # Create MyMotifs instance
    motifs = MyMotifs(sequences)

    # Display results using Streamlit
    st.write("Counts matrix:")
    st.code(motifs.counts)

    st.write("PWM:")
    st.code(motifs.pwm)

    st.write("Sequence alphabet:", motifs.alphabet)

    st.write("Consensus sequence:", motifs.consensus())

    st.write("Masked Consensus sequence:", motifs.masked_consensus())

    st.write("Probability for sequence 'AAACCT':", motifs.probability_sequence("ATGC"))

    st.write("Probability for sequence 'ATACAGAAACCT':", motifs.probability_sequence("ATACAGAAACCT"))


def page1():
    st.title("Deterministic Motif Finding")

    # Input DNA sequences
    sequences = []
    for i in range(1, 5):  # You can adjust the range based on the number of sequences
        sequence = st.text_input(f"Enter Sequence {i}:", key=f"seq_input_{i}")
        sequences.append(sequence.upper())  # Convert to uppercase for consistency

    size = st.slider("Select Motif Size:", min_value=1, max_value=10, value=5)

    # Create an instance of DeterministicMotifFinding
    motif_finder = DeterministicMotifFinding(size, sequences)

    # Perform exhaustive search to find the motif
    motif_indices = motif_finder.exhaustive_search()

    # Display results using Streamlit
    st.subheader("Motif Indices in Each Sequence:")
    st.json({f"Sequence {i + 1}": indices for i, indices in enumerate(motif_indices)})

    st.subheader("Motif Scores:")
    motif_score = motif_finder.score(motif_indices)
    st.write(f"Motif Score: {motif_score}")

    st.subheader("Motifs in Each Sequence:")
    for i in range(len(sequences)):
        k = motif_indices[i]
        motif_sequence = sequences[i][k:k+size]
        st.write(f"Sequence {i + 1}: {motif_sequence}")

def page2():
    st.title("Promoter Region Prediction")
    st.write("An application of DNA BERT2")

    # Sidebar with information
    st.sidebar.header("About")
    st.sidebar.write("This app uses DNA BERT2 to predict the presence of promoter sequence in a given DNA sequence.")

    # User input
    user_input = st.text_area("Enter a DNA sequence:", height=150)

    # Predict when the user provides input
    if st.button("Classify Sequence"):
        if user_input:
            # Call the pred function for prediction
            predicted_class, confidence = pred(user_input)

            # Display the result
            st.subheader("Prediction Result")
            if predicted_class == 1:
                st.success("Promoter Region detected!")
            else:
                st.info("No Promoter Region found.")

            # Display progress bars with percentages
            st.subheader("Class Distribution")
            st.write("1 - Promoter Region found")
            st.progress(confidence)
            st.text(f"{confidence * 100:.2f}%")
            
            st.write("0 - Promoter Region not found")
            st.progress(1 - confidence)
            st.text(f"{(1 - confidence) * 100:.2f}%")

        else:
            st.warning("Please enter a DNA sequence for classification.")


def pred(sequence):
    encoded_input = tokenizer(sequence, return_tensors='pt')
    
    # Pass the encoded input through the model
    with torch.no_grad():
        outputs = model(input_ids=encoded_input['input_ids'], attention_mask=encoded_input['attention_mask'])
        logits = outputs[0]
        predicted_class = logits.argmax(-1).item()
        confidence = logits.softmax(dim=-1)[0, 1].item()

    return predicted_class, confidence

if __name__ == "__main__":
    main()