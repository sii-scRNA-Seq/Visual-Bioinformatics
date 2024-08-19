import { BehaviorSubject, Observable } from 'rxjs';
import { Injectable } from '@angular/core';

import { Block } from './block.interface';
import { DatasetInfo } from './dataset-info';
import { DatasetInfoService } from './dataset-info.service';

@Injectable({
  providedIn: 'root'
})
export class CurrentDatasetService {

  private readonly currentDataset$: BehaviorSubject<DatasetInfo> = new BehaviorSubject<DatasetInfo> ({
    key: '',
    title: '',
    samples: [],
    integration_obs: []
  });
  readonly currentDataset: Observable<DatasetInfo> = this.currentDataset$.asObservable();

  private datasetInfo: DatasetInfo[] = [];

  constructor(private datasetInfoService: DatasetInfoService) {
    this.datasetInfoService.datasetInfo.subscribe(
      (res) => {
        this.datasetInfo = res;
        this.currentDataset$.next(this.datasetInfo[0]);
      },
    );
  }

  onMatSelectValueChanges(block: Block): void {
    if (block.blockId == 'loaddata') {
      const currentDataset = this.datasetInfo.find(dataset => dataset.key == block.parameters[0].value) || this.currentDataset$.getValue();
      this.currentDataset$.next(currentDataset);
    }
  }
}
