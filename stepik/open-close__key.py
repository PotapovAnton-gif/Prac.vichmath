a, b = [int(i) for i in input().split()]
print(a, b)

prod = a*b
func_e = (a - 1)*(b - 1)

for i in range(2, func_e):
    
    if func_e % i != 0:
        d = i
        break

for i in range(2, func_e):

    if  (d*i) % func_e == 1:
        e = i
        break

s_close = (d, prod)
s_open = (e, prod)
print(s_close, s_open, sep='\n',)