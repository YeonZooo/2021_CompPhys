
def fun1():
    a = float(input("a = "))
    b = float(input("b = "))
    return (a**2+b**2)**0.5

def fun2():
    num1 = int(input("num1: ")) 
    num2 = int(input("num2: ")) 
    if num1<num2: min_num = num1
    else: min_num = num2
    for i in range(min_num, 0, -1):
        if num1%i==0 and num2%i==0:
            tmp = i 
            break
    return tmp * num1/tmp * num2/tmp
    
def fun3():
    num1 = int(input("num1: ")) 
    num2 = int(input("num2: ")) 
    if num1<num2: min_num = num1
    else: min_num = num2
    for i in range(min_num, 0, -1):
        if num1%i==0 and num2%i==0:
            result = i 
            break
    return result

def fun4():
    num = str(input("num: "))
    result = 0
    for i in range(len(num)):
        result += int(num[i])
    return result

def fun5():
    myList = [3,5,6,42,65,88,102,201,66]
    isEven = []
    isOdd = []
    for i in myList:
        if i%2==0: isEven.append(i)
        else: isOdd.append(i)
    return {"Even":len(isEven),"Odd":len(isOdd)}

def fun6(m1,m2,m3):
    g = 9.8
    a1 = g*(m1-(m2+m3))/(m1+m2+m3)
    a2 = g*(m2-m3)/(m2+m3) - a1
    a3 = g*(m3-m2)/(m2+m3) - a1
    return {"m1": m1, "a1": round(a1,2),
            "m2": m2, "a2": round(a2,2),
            "m3": m3, "a3": round(a3,2)}

print("1:")
print("c =",fun1())
print("2:")
print(fun2())
print("3:")
print(fun3())
print("4:")
print(fun4())
print("5:")
print(fun5())
print("6-(a):\n", fun6(1,3,2))
print("6-(b):\n", fun6(1,2,2))
print("6-(c):\n", fun6(1,1,1))

