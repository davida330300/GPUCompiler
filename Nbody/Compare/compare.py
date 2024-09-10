import math

class Body:
    def __init__(self, x, y, z, vx, vy, vz):
        self.x = x
        self.y = y
        self.z = z
        self.vx = vx
        self.vy = vy
        self.vz = vz

def euclidean_distance(x1, y1, z1, x2, y2, z2):
    return math.sqrt((x2 - x1) ** 2 + (y2 - y1) ** 2 + (z2 - z1) ** 2)

def velocity_deviation(vx1, vy1, vz1, vx2, vy2, vz2):
    return math.sqrt((vx2 - vx1) ** 2 + (vy2 - vy1) ** 2 + (vz2 - vz1) ** 2)

def compare_files(file1, file2, n_bodies):
    total_pos_deviation = 0.0
    total_vel_deviation = 0.0
    pos_variance = 0.0
    vel_variance = 0.0
    differing_pairs = 0
    
    with open(file1, 'r') as f1, open(file2, 'r') as f2:
        for _ in range(n_bodies):
            line1 = f1.readline().strip()
            line2 = f2.readline().strip()
            
            values1 = list(map(float, line1.split(',')))
            values2 = list(map(float, line2.split(',')))
            
            b1 = Body(*values1)
            b2 = Body(*values2)

            pos_deviation = euclidean_distance(b1.x, b1.y, b1.z, b2.x, b2.y, b2.z)
            vel_deviation = velocity_deviation(b1.vx, b1.vy, b1.vz, b2.vx, b2.vy, b2.vz)
            
            if pos_deviation > 1e-5 or vel_deviation > 1e-5:  # Tolerance for floating-point precision
                differing_pairs += 1
                total_pos_deviation += pos_deviation
                total_vel_deviation += vel_deviation
                pos_variance += pos_deviation ** 2
                vel_variance += vel_deviation ** 2
    
    if differing_pairs == 0:
        print("The files are identical.")
    else:
        avg_pos_deviation = total_pos_deviation / differing_pairs
        avg_vel_deviation = total_vel_deviation / differing_pairs
        pos_stdev = math.sqrt(pos_variance / differing_pairs)
        vel_stdev = math.sqrt(vel_variance / differing_pairs)
        
        print(f"The files differ in {differing_pairs} pairs.")
        print(f"Average position deviation: {avg_pos_deviation:.6f}")
        print(f"Standard deviation of position deviation: {pos_stdev:.6f}")
        print(f"Average velocity deviation: {avg_vel_deviation:.6f}")
        print(f"Standard deviation of velocity deviation: {vel_stdev:.6f}")

if __name__ == "__main__":
    file1 = "gpu_end.csv"
    file2 = "serial_end.csv"
    n_bodies = 10000  # Adjust according to your data
    compare_files(file1, file2, n_bodies)
