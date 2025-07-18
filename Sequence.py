from print_string_sequence import print_string_sequence

class Sequence():
    def __init__(self, data_path):
        self.Data = ""
        with open(data_path, "r") as file1:
            for line in file1.readlines():
                if line.startswith('>'):
                    self.Header = line[:-2]
                else:
                    self.Data = self.Data + line[:-2]
        self.DataSize = len(self.Data)

    def getItr(self, itr):
        return self.Data[itr]
    
    def getSize(self):
        return self.DataSize
    
    def find(self, ch):
        return self.Data.rfind(ch)
        
    def print(self):
        print_string_sequence(self.Data)