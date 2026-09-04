"""TensorFlow/Keras implementation of the parity benchmark."""

from __future__ import annotations

import random
from pathlib import Path

import numpy as np

from .common import (
    BenchmarkTimer,
    build_vocabulary,
    classification_metrics,
    error_analysis,
    vectorize,
)
from .dataset import LABELS, DatasetBundle


def run_tensorflow(bundle: DatasetBundle, output_dir: Path, seed: int = 2026, epochs: int = 40, patience: int = 5) -> dict[str, object]:
    try:
        import tensorflow as tf
    except ImportError as exc:
        raise RuntimeError("TensorFlow dependency blocker: install tensorflow to run the TensorFlow benchmark") from exc

    random.seed(seed)
    np.random.seed(seed)
    tf.keras.utils.set_random_seed(seed)
    try:
        tf.config.experimental.enable_op_determinism()
    except (AttributeError, RuntimeError):
        pass
    vocabulary = build_vocabulary(bundle.train)
    x_train, y_train = vectorize(bundle.train, vocabulary)
    x_val, y_val = vectorize(bundle.validation, vocabulary)
    x_test, y_test = vectorize(bundle.test, vocabulary)
    train_ds = tf.data.Dataset.from_tensor_slices((x_train, y_train)).shuffle(len(x_train), seed=seed, reshuffle_each_iteration=False).batch(16)
    model = tf.keras.Sequential([
        tf.keras.layers.Input(shape=(len(vocabulary),)),
        tf.keras.layers.Dense(24, activation="relu"),
        tf.keras.layers.Dense(len(LABELS)),
    ])
    model.compile(optimizer=tf.keras.optimizers.Adam(0.02), loss=tf.keras.losses.SparseCategoricalCrossentropy(from_logits=True))
    output_dir.mkdir(parents=True, exist_ok=True)
    checkpoint = output_dir / "tensorflow_checkpoint.weights.h5"
    callbacks = [
        tf.keras.callbacks.EarlyStopping(monitor="val_loss", patience=patience, restore_best_weights=True),
        tf.keras.callbacks.ModelCheckpoint(checkpoint, monitor="val_loss", save_best_only=True, save_weights_only=True),
    ]
    with BenchmarkTimer() as timer:
        history_obj = model.fit(train_ds, validation_data=(x_val, y_val), epochs=epochs, callbacks=callbacks, verbose=0)
    predicted = np.argmax(model.predict(x_test, verbose=0), axis=1)
    history = [
        {"epoch": index + 1, "validation_loss": float(loss)}
        for index, loss in enumerate(history_obj.history["val_loss"])
    ]
    metrics = classification_metrics(y_test, predicted)
    return {
        "framework": "tensorflow", "framework_version": tf.__version__, "seed": seed,
        "parameter_count": model.count_params(), "epochs_completed": len(history),
        "best_validation_loss": min(history_obj.history["val_loss"]), "runtime_seconds": timer.seconds,
        "python_peak_memory_bytes": timer.peak_bytes, "checkpoint": str(checkpoint),
        "metrics": metrics, "errors": error_analysis(bundle.test, y_test, predicted), "history": history,
    }
