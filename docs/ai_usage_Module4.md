- AI Use Log
- Tool/model & version:
- What I asked for:
- Snippet of prompt(s):
- What I changed before committing:
- How I verified correctness (tests, sample data):

- Tool/model & version: ChatGPT 5
- What I asked for: Pseudocode breakdown
- Snippet of prompt(s): Breakdown this pseudocode for me: PatternCount(Text, Pattern) count ← 0 for i ← 0 to |Text| − |Pattern| if Text(i, |Pattern|) = Pattern count ← count + 1 return count
- What I changed before committing: N/A
- How I verified correctness (tests, sample data): N/A

- Tool/model & version: ChatGPT 5
- What I asked for: Python code for the above
- Snippet of prompt(s): Please give me the python code to run. 
- What I changed before committing: Added def pattern_count(text, pattern):
    count = 0
    # Loop through each possible starting position
    for i in range(len(text) - len(pattern) + 1):
        # Slice the substring of the same length as pattern
        if text[i:i+len(pattern)] == pattern:
            count += 1
    return count
- How I verified correctness (tests, sample data): N/A yet

- Tool/model & version: ChatGPT 5
- What I asked for: Python code to use it
- Snippet of prompt(s): Please give me the python code for an example use. 
- What I changed before committing: Added text, pattern, print.
- How I verified correctness (tests, sample data): Did not print anything.

- Tool/model & version: Chat GPT 4
- What I asked for: Why my code didn't print
- Snippet of prompt(s): Why didn't this code print? 
- What I changed before committing: def pattern_count(text, pattern): count = 0 # Loop through each possible starting position for i in range(len(text) - len(pattern) + 1): # Slice the substring of the same length as pattern if text[i:i+len(pattern)] == pattern: count += 1 return count text = GCGCG pattern = GCG print(pattern_count(text, pattern))
- How I verified correctness (tests, sample data): Code printed this time. Verified against Rosalind data.

- Tool/model & version: Chat GPT 4
- What I asked for: Code check
- Snippet of prompt(s): Sometimes this code works, sometimes it does not. Can you analyze why?
- What I changed before committing: Added # Read DNA string and clean it
with open(filename, 'r') as file:
    dna_text = "".join(file.read().split()).upper()  # remove whitespace/newlines and capitalize
Set k for each one. 
- How I verified correctness (tests, sample data): Checked against Rosalind debugging sets.

- Tool/model & version: ChatGPT 4
- What I asked for: Problem in plain English.
- Snippet of prompt(s): Please restate the following Rosalind problem in plain english, identify inputs and outputs, and describe an algorithm to solve it, but do not give me code yet.
Given: A string Genome, and integers k, L, and t.
Return: All distinct k-mers forming (L, t)-clumps in Genome.
- What I changed before committing: N/A
- How I verified correctness (tests, sample data): N/A

- Tool/model & version: ChatGPT 4
- What I asked for: Psuedocode
- Snippet of prompt(s): Please put the above problem into pseudocode. 
- What I changed before committing: N/A
- How I verified correctness (tests, sample data): N/A

- Tool/model & version: ChatGPT 4
- What I asked for: Python code
- Snippet of prompt(s): Please turn the above pseduocode into python for Colab.
- What I changed before committing: Added code set. 
- How I verified correctness (tests, sample data): Ran against Rosalind debugging sets.

- Tool/model & version: ChatGPT 4
- What I asked for: Problem in plain English.
- Snippet of prompt(s): Please restate the following Rosalind problem in plain english, identify inputs and outputs, and describe an algorithm to solve it, but do not give me code yet. Given: A DNA string Genome. Return: All integer(s) i minimizing Skew(Prefixi (Text)) over all values of i (from 0 to |Genome|).
- What I changed before committing: N/A
- How I verified correctness (tests, sample data): N/A

- Tool/model & version: ChatGPT 4
- What I asked for: Psuedocode
- Snippet of prompt(s): Please put the above problem into pseudocode. 
- What I changed before committing: N/A
- How I verified correctness (tests, sample data): N/A

- Tool/model & version: ChatGPT 4
- What I asked for: Python code
- Snippet of prompt(s): Please turn the above pseduocode into python for Colab.
- What I changed before committing: Added code set. 
- How I verified correctness (tests, sample data): Ran against Rosalind debugging sets.

- Tool/model & version: ChatGPT 4
- What I asked for: Understanding check
- Snippet of prompt(s): Ask me 2 questions to check my understanding without code. 
- What I changed before committing: N/A 
- How I verified correctness (tests, sample data): N/A

- Tool/model & version: ChatGPT 4
- What I asked for: Problem in plain English.
- Snippet of prompt(s): Please restate the following Rosalind problem in plain english, identify inputs and outputs, and describe an algorithm to solve it, but do not give me code yet. Given: Two DNA strings. Return: An integer value representing the Hamming distance.
- What I changed before committing: N/A
- How I verified correctness (tests, sample data): N/A

- Tool/model & version: ChatGPT 4
- What I asked for: Psuedocode
- Snippet of prompt(s): Please put the above problem into pseudocode. 
- What I changed before committing: N/A
- How I verified correctness (tests, sample data): N/A

- Tool/model & version: ChatGPT 4
- What I asked for: Python code
- Snippet of prompt(s): Please turn the above pseduocode into python for Colab.
- What I changed before committing: Added code set. 
- How I verified correctness (tests, sample data): Ran against Rosalind debugging sets.

- Tool/model & version: ChatGPT 4
- What I asked for: Understanding check
- Snippet of prompt(s): Ask me 2 questions to check my understanding without code. 
- What I changed before committing: N/A 
- How I verified correctness (tests, sample data): N/A

- Tool/model & version: ChatGPT 4
- What I asked for: Problem in plain English.
- Snippet of prompt(s): Please restate the following Rosalind problem in plain english, identify inputs and outputs, and describe an algorithm to solve it, but do not give me code yet. Given: A DNA string Pattern and an integer d. Return: The collection of strings Neighbors(Pattern, d).
- What I changed before committing: N/A
- How I verified correctness (tests, sample data): N/A

- Tool/model & version: ChatGPT 4
- What I asked for: Psuedocode
- Snippet of prompt(s): Please put the above problem into pseudocode. 
- What I changed before committing: N/A
- How I verified correctness (tests, sample data): N/A

- Tool/model & version: ChatGPT 4
- What I asked for: Python code
- Snippet of prompt(s): Please turn the above pseduocode into python for Colab.
- What I changed before committing: Added code set. 
- How I verified correctness (tests, sample data): Ran against Rosalind sample set.

- Tool/model & version: ChatGPT 4
- What I asked for: Print layout.
- Snippet of prompt(s): How do I get this to print in a column?
- What I changed before committing: Added # Sort the neighbors alphabetically
sorted_neighbors = sorted(neighbors)
# Print each neighbor on a new line (column format)
for neighbor in sorted_neighbors:
    print(neighbor)
- How I verified correctness (tests, sample data): Ran against Rosalind debug sets.

- Tool/model & version: ChatGPT 4
- What I asked for: Understanding check
- Snippet of prompt(s): Ask me 2 questions to check my understanding without code. 
- What I changed before committing: N/A 
- How I verified correctness (tests, sample data): N/A

- Tool/model & version: ChatGPT 4
- What I asked for: Problem in plain English.
- Snippet of prompt(s): Please restate the following Rosalind problem in plain english, identify inputs and outputs, and describe an algorithm to solve it, but do not give me code yet. Given: A string Text as well as integers k and d. Return: All most frequent k-mers with up to d mismatches in Text.
- What I changed before committing: N/A
- How I verified correctness (tests, sample data): N/A

- Tool/model & version: ChatGPT 4
- What I asked for: Psuedocode
- Snippet of prompt(s): Please put the above problem into pseudocode. 
- What I changed before committing: N/A
- How I verified correctness (tests, sample data): N/A

- Tool/model & version: ChatGPT 4
- What I asked for: Python code
- Snippet of prompt(s): Please turn the above pseduocode into python for Colab.
- What I changed before committing: Added code set. 
- How I verified correctness (tests, sample data): Did not run correctly. 

- Tool/model & version: ChatGPT 4
- What I asked for: Code update
- Snippet of prompt(s): Please update the code since Rosalind has k and d on the same line. 
- What I changed before committing: Added # Rosalind input format:
# Line 1: Text
# Line 2: k d (space-separated)
Text = lines[0].strip()
k, d = map(int, lines[1].strip().split())  # Read both integers from the same line
- How I verified correctness (tests, sample data): Ran against Rosalind debug sets.

- Tool/model & version: ChatGPT 4
- What I asked for: Understanding check
- Snippet of prompt(s): Ask me 2 questions to check my understanding without code. 
- What I changed before committing: N/A 
- How I verified correctness (tests, sample data): N/A

- Tool/model & version: ChatGPT 4
- What I asked for: Problem in plain English.
- Snippet of prompt(s): Please restate the following Rosalind problem in plain english, identify inputs and outputs, and describe an algorithm to solve it, but do not give me code yet. Given: A DNA string Text as well as integers k and d. Return: All k-mers Pattern maximizing the sum Countd(Text, Pattern) + Countd(Text, Pattern) over all possible k-mers.
- What I changed before committing: N/A
- How I verified correctness (tests, sample data): N/A

- Tool/model & version: ChatGPT 4
- What I asked for: Psuedocode
- Snippet of prompt(s): Please put the above problem into pseudocode. 
- What I changed before committing: N/A
- How I verified correctness (tests, sample data): N/A

- Tool/model & version: ChatGPT 4
- What I asked for: Python code
- Snippet of prompt(s): Please turn the above pseduocode into python for Colab.
- What I changed before committing: Added code set. 
- How I verified correctness (tests, sample data): Did not run correctly.

- Tool/model & version: ChatGPT 4
- What I asked for: Understanding check
- Snippet of prompt(s): Ask me 2 questions to check my understanding without code. 
- What I changed before committing: N/A 
- How I verified correctness (tests, sample data): N/A
