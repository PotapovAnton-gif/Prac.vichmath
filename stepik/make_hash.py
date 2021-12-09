import hashlib
  
# initializing string
st = ["Blue", "Grass", "Sun", "HelloHome", "Привет", "Sky", "привет"]

  
for i in st:
    result = hashlib.sha512(i.encode())
    
    # printing the equivalent hexadecimal value.
    print("The hexadecimal equivalent of SHA256 is : {}".format(i))
    print(result.hexdigest())
    
    print ("\r")