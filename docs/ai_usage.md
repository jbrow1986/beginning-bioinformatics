# AI Use Log
- Tool/model & version:
- What I asked for:
- Snippet of prompt(s):
- What I changed before committing:
- How I verified correctness (tests, sample data):

- Tool/model & version: ChatGPT 5
- What I asked for: Help defining functions and arguments.
- Snippet of prompt(s): In coding python, parentheses are used to provide arguments to functions. Can you break that down for me?
- What I changed before committing: N/A
- How I verified correctness (tests, sample data): N/A

- Tool/model & version: ChatGPT 5
- What I asked for: Counter variable explanation
- Snippet of prompt(s): In python coding, what is a counter variable?
- What I changed before committing: N/A
- How I verified correctness (tests, sample data): N/A

- Tool/model & version: ChatGPT 5
- What I asked for: Code breakdown
- Snippet of prompt(s): Can you break this code down for me? with open("practice.txt", "r") as data: line_number = 0 for line in data: if line_number % 2 == 0: # Check if the line number is even (starting from 0) print(line.rstrip()) line_number += 1
- What I changed before committing: Nothing
- How I verified correctness (tests, sample data): Checked against Rosalind practice data.

- Tool/model & version: ChatGPT 5
- What I asked for: Code check
- Snippet of prompt(s): I'm trying to code to read every other line of a file. What's wrong with my code? with open("practice.txt","r") as data: line_number = 1 for line in data: if line_number % 2 == 0: print (list_of_lines) line_number += 1
- What I changed before committing: Changed print (list_of_lines) to print(line.rstrip()). Changed indent on line_number.
- How I verified correctness (tests, sample data): Tested against Rosalind data sample.

- Tool/model & version: ChatGPT 5
- What I asked for: Code check.
- Snippet of prompt(s): I'm trying to get my code to read every other line using the join function. Can you tell me what's wrong with my code? with open("practice.txt","r") as data: list_of_lines=data.readlines() line_number = 1 if line_number % 2 == 0: print ("".join(list_of_lines))
- What I changed before committing: Changed back to for loop instead of join slicing.
- How I verified correctness (tests, sample data): Checked against Rosalind dataset.

- Tool/model & version: ChatGPT 5
- What I asked for: Explain the difference between 2 sets of code.
- Snippet of prompt(s): What is the difference between these 2 sets of code? with open("practice.txt","r") as data: line_number = 1 for line in data: if line_number % 2 == 0: print(line.rstrip()) line_number += 1 with open("practice.txt","r") as data: list_of_lines=data.readlines() line_number = 1 every_other = [] for line in list_of_lines: if line_number % 2 == 0: every_other.append(line) line_number += 1 print ("".join(every_other))
- What I changed before committing: N/A
- How I verified correctness (tests, sample data): N/A

- Tool/model & version: ChatGPT 5
- What I asked for: Code check
- Snippet of prompt(s): What's wrong with this code? with open("/content/rosalinddata.txt","r") as rosalinddata: line_number = 1 for line in rosalinddata: if line_number % 2 == 0: print(line.rstrip()) line_number += 1
- What I changed before committing: Changed indentation on the for loop.
- How I verified correctness (tests, sample data): Checked against Rosalind data.

- Tool/model & version: ChatGPT 5
- What I asked for: How do I break a dictionary down to print line by line?
- Snippet of prompt(s): How do I break a dictionary down to print line by line?
- What I changed before committing: changed lines for print key and print value.
- How I verified correctness (tests, sample data): Manually checked it printed correctly.

- Tool/model & version: ChatGPT5
- What I asked for: Coding read and loop at the same time.
- Snippet of prompt(s): If I've already coded for a file to be read, can I then loop over it?
- What I changed before committing: Added lines = data.readlines()   # puts all lines in a list
    for line in lines:
        print(line.rstrip())
- How I verified correctness (tests, sample data): Manually checked it printed correctly.

- Tool/model & version: ChatGPT 5
- What I asked for: Finding words one at a time.
- Snippet of prompt(s): How do I loop over a file to get words one at a time?
- What I changed before committing: Added with open("practice.txt", "r") as data:
    for line in data:              
        words = line.split()       
        for word in words:         
            print(word)
- How I verified correctness (tests, sample data): Manually checked it printed correctly.

- Tool/model & version: Gemini 1.5 Flash
- What I asked for: Explain an attribute error
- Snippet of prompt(s): Please explain this error: Please explain this error: AttributeError: '_io.TextIOWrapper' object has no attribute 'count'
- What I changed before committing: Added DNA = to practice 6 line. 
- How I verified correctness (tests, sample data): Tested against Rosalind data.

- Tool/model & version: Gemini 1.5 Flash
- What I asked for: Error check
- Snippet of prompt(s): What does this error mean: AttributeError: 'dict' object has no attribute 'append'
- What I changed before committing: Nothing.
- How I verified correctness (tests, sample data): Did not change due to Gemini's answer.

- Tool/model & version: ChatGPT 5
- What I asked for: Code check.
- Snippet of prompt(s): Help me with this code. I'm trying to count unique words in file. with open("/content/rosalinddata5.txt","r") as rosalinddata5: dictionary = {} for line in rosalinddata5: words = line.split() for word in words: if word not in dictionary: dictionary.append(word) for key, value in dictionary.items(): print (count(key, value))
- What I changed before committing: Added: if word not in word_counts:
                word_counts[word] = 1    
            else:
                word_counts[word] += 1
  for word, count in word_counts.items():
for word, count in word_counts.items():
    print(word, count)
- How I verified correctness (tests, sample data): Compared to previous question. Matched that code, not the answer I was looking for.

- Tool/model & version: ChatGPT 5
- What I asked for: Counting unique words.
- Snippet of prompt(s): This is not quite what I'm looking for. I need the number of unique words in the entire dictionary. So it should count each word only once and then give me a summary.
- What I changed before committing: Deleted if word not in dictionary:
        dictionary[word] = 1
      else:
        dictionary[word] += 1
for key, value in dictionary.items():
  print (key, value).
  Added: word_count.add(word)
count = 0
for word in word_count:
  count += 1
print (count)
- How I verified correctness (tests, sample data): Tested with a document of 5 unique words. 
