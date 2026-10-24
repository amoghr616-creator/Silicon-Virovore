def solve(A, B, cA, cB, fA):
    init = B // cB * cA + A
    if init >= fA:
        return 0
    nA0 = fA - 1 - init
    y = cB - 1 - B % cB  # nB0
    if cA >= cB:
        y += nA0
    else:
        y += nA0 // cA * cB + nA0 % cA
    return y + 1

def main():
    T = int(input())
    for _ in range(T):
        A, B, cA, cB, fA = map(int, input().split())
        print(solve(A, B, cA, cB, fA))

if __name__ == "__main__":
    main()