import os
import cv2
import numpy as np

def overlay_colored_mask_on_blue(image_dir, mask_dir, output_dir="pseudocolored", alpha=1):
    os.makedirs(output_dir, exist_ok=True)

    # Build mask lookup dictionary from filenames without "_Clusters"
    mask_files = [f for f in os.listdir(mask_dir) if f.lower().endswith(('.tif', '.tiff', '.png', '.jpg', '.jpeg'))]
    mask_dict = {os.path.splitext(f.replace('_Clusters', ''))[0]: f for f in mask_files}

    # Loop through each image file
    image_files = [f for f in os.listdir(image_dir) if f.lower().endswith(('.tif', '.tiff', '.png', '.jpg', '.jpeg'))]

    for img_file in image_files:
        base_name = os.path.splitext(img_file)[0]

        if base_name not in mask_dict:
            print(f"No matching mask for {img_file}, skipping.")
            continue

        image_path = os.path.join(image_dir, img_file)
        mask_path = os.path.join(mask_dir, mask_dict[base_name])

        # Load RGB image
        image = cv2.imread(image_path)
        if image is None:
            print(f"Failed to read image: {img_file}")
            continue

        # Keep only the blue channel (others set to 0)
        blue_only = np.zeros_like(image)
        blue_only[:, :, 0] = image[:, :, 0]  # keep blue, zero out green and red

        # Load the colored mask (keep all 3 channels)
        mask = cv2.imread(mask_path, cv2.IMREAD_COLOR)
        if mask is None:
            print(f"Failed to read mask: {mask_path}")
            continue

        # Create a binary mask from non-black pixels in the mask
        mask_binary = np.any(mask > 0, axis=-1)  # shape: (H, W)

        # Prepare for blending
        blended = blue_only.copy()

        # Only blend where the mask has content
        for c in range(3):
            blended[:, :, c] = np.where(
                mask_binary,
                (1 - alpha) * blue_only[:, :, c] + alpha * mask[:, :, c],
                blue_only[:, :, c]
            ).astype(np.uint8)

        # Save result
        output_path = os.path.join(output_dir, img_file)
        cv2.imwrite(output_path, blended)
        print(f"Saved: {output_path}")


# Example usage
overlay_colored_mask_on_blue(
    image_dir="/Users/alexlawson/Masters-Data-Final/Representative-images/PVN-02/original/renamed_output",
    mask_dir="/Users/alexlawson/Masters-Data-Final/Morphoglia/PVN_all/Morphoglia/MorphoGlia/ID_clusters", 
    output_dir="/Users/alexlawson/Masters-Data-Final/Representative-images/Morphoglia/pseudo-v3"
)
