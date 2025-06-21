<img src="docs/logo.png" width="300"/>

# Why We’re Still Here: A Simulation of Human Evolution, War, and the Power of Getting Along

##### by Francesco Giuseppe Gillio
###### Department of Control and Computer Engineering
###### Politecnico di Torino

## Reference Paper
- [Why We’re Still Here: A Simulation of Human Evolution, War, and the Power of Getting Along](https://github.com/sencosx/natural-selection-simulator/blob/main/docs/natural-selection-simulator-ACM.pdf)

---

## Table of Contents

- [What’s the Problem? (Abstract)](#whats-the-problem-abstract)
- [The Tools of the Trade (Algorithms)](#the-tools-of-the-trade-algorithms)

---

## What’s the Problem? (Abstract)

This study builds a discrete-event simulator with a clear purpose: to recreate the messy, competitive reality of evolution when multiple Homo species shared the same terrain and fought—knowingly or not—for survival. The simulation zooms in on the early Holocene, about 11,650 years ago, to test what might have happened when Homo sapiens, Homo neanderthalensis, and Homo erectus overlapped in space and time.
    
Instead of guessing or speculating, this work runs the tape forward using high-performance computational models. It doesn't just ask who survived, but why. The simulation tests the role of social structures, cognitive differences, and competition for limited resources. The results are sharp: cooperation beats isolation, and strategy matters more than brute strength.
    
What emerges is not a romantic tale of survival, but a data-driven look at extinction. Social alliances, flexible thinking, and collective behavior turn out to be evolutionary game-changers. This project doesn’t just add to our understanding of ancient human dynamics—it sets the stage for a new way of studying evolution itself: not as a slow drift, but as a dynamic system where complexity, pressure, and choice collide.

---

## The Tools of the Trade (Algorithms)

Play time starts here.

1. Fire this up in your notebook:

```python
!git clone https://github.com/sencosx/natural-selection-simulator.git
```

2. Run the simulator and see what happens:

```python
!python /content/natural-selection-simulator/src/simulator.py \
        --selections 10000 \
        --timeframe 45 \
        --rate 0.15 \
        --alpha 0.25 \
        --improve 0.5 \
        --lifetime 15 \
        --directory output
```

Boom. **Game over**.

---

> *“If your code always makes sense to you, you're either an idiot — or not thinking hard enough.”*

That’s what this project is about. Now go dig in.

— *F.G.G. (still writing code, as always)*
