from classify_protein import score_sequence
import numpy as np



def test1():
  query = "A"
  transition_probs = np.log(np.array([
    # S,  I_0,    M_1,     D_1
    [ 0,   1.0,    0.0,     0.0],  # S
    [ 0,   0.0,    1.0,     0.0],  # I_0
    [ 0,   0.0,    0.0,     1.0],  # M_1
    [ 0,   0.0,    0.0,     0.0],  # D_1
  ]))

  m_emission_probs = np.log(np.array([
    # A    R    N    D    C
    [ 0.0, 0.0, 0.0, 0.0, 0.0],  # M_0 (does not exist)
    [ 1.0, 0.0, 0.0, 0.0, 0.0],  # M_1
  ]))

  i_emission_probs = np.log(np.array([
    # A    R    N    D    C
    [ 1.0, 0.0, 0.0, 0.0, 0.0],  # I_0
    [ 1.0, 0.0, 0.0, 0.0, 0.0],  # I_1
  ]))

  max_length = 1

  # query "A" matches exactly with state M_1
  path, prob = score_sequence(query, transition_probs, m_emission_probs, i_emission_probs, max_length)
  path = "M"
  prob = 0  

  try:
      assert path == "M", f"Expected path to be 'M', but got {path}"
      
      assert prob == 0, f"Expected score to be 0, but got {prob}"
      
      print("Test 1 Correct!")
      
  except AssertionError as e:
      print(f"AssertionError: {e}")

test1()


import numpy as np

def test2():
    query = "ARN" 
    transition_probs = np.array([
        # S,   I_0,   M_1,   D_1,   I_1,   M_2,   D_2,   I_2,   M_3,   D_3,   I_3
        [ 0,   0.0,   1.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0],  # S
        [ 0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0],  # I_0
        [ 0,   0.0,   0.0,   0.0,   0.0,   0.0,   1.0,   0.0,   0.0,   0.0,   0.0],  # M_1
        [ 0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0],  # D_1
        [ 0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0],  # I_1
        [ 0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0],  # M_2
        [ 0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   1.0,   0.0],  # D_2
        [ 0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0],  # I_2
        [ 0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0],  # M_3
        [ 0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0],  # D_3
        [ 0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0],  # I_3
        
    ])


    m_emission_probs = np.array([
        #   A,   R,   N,   D,   C
        [ 0.0, 0.0, 0.0, 0.0, 0.0],  # M_0 (does not exist)
        [ 1.0, 0.0, 0.0, 0.0, 0.0],  # M_1
        [ 0.0, 0.0, 0.0, 0.0, 0.0],  # M_2
        [ 0.0, 0.0, 0.0, 0.0, 0.0],  # M_3
    ])

    i_emission_probs = np.array([
        #   A,   R,   N,   D,   C
        [ 0.0, 0.0, 0.0, 0.0, 0.0],  
        [ 0.0, 0, 0, 0.0, 0.0],  # I_1 emits both R & N
    ])

    max_length = 1

    # query is "ARN"
    path, prob = score_sequence(query, transition_probs, m_emission_probs, i_emission_probs, max_length)
    
    # "A" is emitted by M_1 -> D2 -> D3
    expected_path = "MDD"
    expected_prob = np.log(1.0) * 3 

    try:
        assert path == expected_path, f"Expected path to be {expected_path}, but got {path}"
        assert prob == expected_prob, f"Expected score to be {expected_prob}, but got {prob}"
        print("Test 2 Correct!")
    except AssertionError as e:
        print(f"AssertionError: {e}")

def test3():
    # insertions
    print("test 3 insetionts")
    query = "C"
    
    transition_probs = np.array([
        # S,    I_0,    M_1,    D_1,    I_1,    M_2,    D_2,    I_2,    M_3,    D_3,    I_3
        [ 0,    0.0,    1.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0],  # S
        [ 0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0],  # I_0
        [ 0,    0.0,    0.0,    0.0,    1.0,    0.0,    0.0,    1.0,    0.0,    0.0,    0.0],  # M_1
        [ 0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0],  # D_1
        [ 0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    1.0,    0.0,    0.0,    0.0],  # I_1
        [ 0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0],  # M_2
        [ 0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0],  # D_2
        [ 0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    1.0],  # I_2
        [ 0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0],  # M_3
        [ 0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0],  # D_3
        [ 0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    1.0],  # I_3
    ])

    m_emission_probs = np.array([
        # A    R    N    D    C
        [ 0.0, 0.0, 0.0, 0.0, 0.0],  # M_0 (does not exist)
        [ 0.0, 0.0, 0.0, 0.0, 1.0],  # M_1 = C
        [ 0.0, 0.0, 0.0, 0.0, 0.0],  # M_2
        [ 0.0, 0.0, 0.0, 0.0, 0.0],  # M_3
    ])

    i_emission_probs = np.array([
        # A    R    N    D    C
        [ 0.0, 0.0, 0.0, 0.0, 0.0],  # I_0
        [ 0.0, 0.0, 0.0, 0.0, 0.0],  # I_1
        [ 0.0, 0.0, 1.0, 0.0, 0.0],  # I_2
        [ 0.0, 0.0, 0.0, 1.0, 0.0],  # I_3
    ])

    max_length = 3

    # C is emitted by M_1 and then insertions I_2 -> I_2
    expected_path = "MII"
    expected_prob = np.log(1.0) + np.log(1.0)  + np.log(1.0)

    path, prob = score_sequence(query, transition_probs, m_emission_probs, i_emission_probs, max_length)

    try:
        assert path == expected_path, f"Expected path to be {expected_path}, but got {path}"
        assert prob == expected_prob, f"Expected score to be {expected_prob}, but got {prob}"
        print("Test 3 Correct!")
    except AssertionError as e:
        print(f"AssertionError: {e}")



def test4():
    query = "CR"
    
    transition_probs = np.array([
        # S,    I_0,    M_1,    D_1,    I_1,    M_2,    D_2,    I_2,    M_3,    D_3,    I_3,    M_4,    D_4,    I_4
        [ 0,    0.0,    1.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0],  # S
        [ 0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0],  # I_0
        [ 0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    1.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0],  # M_1
        [ 0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0],  # D_1
        [ 0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0],  # I_1
        [ 0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0],  # M_2
        [ 0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    1.0,    0.0,    0.0,    0.0,    0.0],  # D_2
        [ 0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    1.0,    0.0,    0.0,    0.0],  # I_2
        [ 0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0],  # M_3
        [ 0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0],  # D_3
        [ 0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    1.0,    0.0,    0.0],  # I_3
        [ 0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0],  # M_4
        [ 0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0],  # D_4
        [ 0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0],  # I_4
    ])

    m_emission_probs = np.array([
        # A    R    N    D    C
        [ 0.0, 0.0, 0.0, 0.0, 0.0],  # M_0 (does not exist)
        [ 0.0, 0.0, 0.0, 0.0, 1.0],  # M_1 = C
        [ 0.0, 0.0, 0.0, 0.0, 0.0],  # M_2
        [ 0.0, 0.0, 0.0, 0.0, 0.0],  # M_3
        [ 0.0, 1.0, 0.0, 0.0, 0.0],  # M_4 = R
    ])

    i_emission_probs = np.array([
        # A    R    N    D    C
        [ 0.0, 0.0, 0.0, 0.0, 0.0],  # I_0
        [ 0.0, 0.0, 0.0, 0.0, 0.0],  # I_1
        [ 0.0, 0.0, 1.0, 0.0, 0.0],  # I_2 = N
        [ 0.0, 0.0, 0.0, 1.0, 0.0],  # I_3 = D
        [ 0.0, 0.0, 0.0, 0.0, 0.0],  # I_4
    ])

    max_length = 4

    expected_path = "MIIM"
    expected_prob = np.log(1.0)

    path, prob = score_sequence(query, transition_probs, m_emission_probs, i_emission_probs, max_length)

    try:
        assert path == expected_path, f"Expected path to be {expected_path}, but got {path}"
        assert prob == expected_prob, f"Expected score to be {expected_prob}, but got {prob}"
        print("Test 4 Correct!")
    except AssertionError as e:
        print(f"AssertionError: {e}")

def main():
    test1() 
    test2() #deletions at end
    test3() # insertions at end 
    test4() #insertions in middle

main()

