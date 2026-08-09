import os

def get_poscar_atom_count(poscar_path):
    """POSCAR 파일에서 총 원자 수를 읽어옵니다."""
    with open(poscar_path, 'r') as f:
        lines = f.readlines()
    
    try:
        # 6번째 줄(index 5) 확인: 알파벳(원소 기호)이 포함되어 있으면 VASP 5 포맷
        line_5 = lines[5].strip()
        if any(c.isalpha() for c in line_5):
            counts = lines[6].split()  # 원자 수는 7번째 줄(index 6)에 있음
        else:
            counts = line_5.split()    # VASP 4 포맷이면 6번째 줄(index 5)이 원자 수
            
        return sum(int(x) for x in counts)
    except Exception as e:
        print(f"POSCAR 파싱 에러 ({poscar_path}): {e}")
        return None

# 처음 12개 원자에 대한 고정 MAGMOM 설정
BASE_MAGMOM = "2*3.8000 2*-3.8000 2*3.8000 2*-3.8000 4*4.4"
BASE_ATOMS = 12

# 현재 디렉토리 내의 모든 하위 폴더 탐색
subdirs = [d for d in os.listdir('.') if os.path.isdir(d)]

for d in subdirs:
    incar_path = os.path.join(d, 'INCAR')
    poscar_path = os.path.join(d, 'POSCAR')

    # INCAR와 POSCAR가 모두 존재하는 폴더만 처리
    if not (os.path.exists(incar_path) and os.path.exists(poscar_path)):
        continue

    total_atoms = get_poscar_atom_count(poscar_path)
    if total_atoms is None:
        continue

    with open(incar_path, 'r') as f:
        incar_lines = f.readlines()

    # 기존 INCAR에서 MAGMOM과 NCORE가 적힌 줄은 모두 필터링하여 제거 (중복 방지)
    new_incar_lines = []
    for line in incar_lines:
        line_upper = line.strip().upper()
        if not (line_upper.startswith('MAGMOM') or line_upper.startswith('NCORE')):
            new_incar_lines.append(line)
    
    # 추가해야 할 원자 수 계산
    diff = total_atoms - BASE_ATOMS
    
    if diff > 0:
        new_magmom_line = f"MAGMOM = {BASE_MAGMOM} {diff}*0.0000\n"
    elif diff == 0:
        new_magmom_line = f"MAGMOM = {BASE_MAGMOM}\n"
    else:
        print(f"[{d}] 경고: 총 원자 수({total_atoms})가 12개 미만이라 건너뜁니다.")
        continue

    # 파일의 맨 끝에 새로운 MAGMOM과 NCORE 추가
    new_incar_lines.append(new_magmom_line)
    new_incar_lines.append("NCORE = 8\n")

    # INCAR 덮어쓰기
    with open(incar_path, 'w') as f:
        f.writelines(new_incar_lines)
        
    print(f"[{d}] 업데이트 완료 (MAGMOM 세팅 및 NCORE = 8 추가)")
