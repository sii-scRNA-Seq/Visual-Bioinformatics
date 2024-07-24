/* eslint-disable @typescript-eslint/no-unused-vars */

import { Injectable } from '@angular/core';
import { BehaviorSubject, Observable } from 'rxjs';

import { Block } from './block.interface';
import { Output } from './output';
import { OutputServiceInterface } from './output.service.interface';

@Injectable({
  providedIn: 'root'
})
export class MockOutputService implements OutputServiceInterface {
  private readonly executingBlocks$: BehaviorSubject<boolean> = new BehaviorSubject<boolean> (false);
  readonly executingBlocks: Observable<boolean> = this.executingBlocks$.asObservable();

  private readonly outputs$: BehaviorSubject<Output[]> = new BehaviorSubject<Output[]> ([]);
  readonly outputs: Observable<Output[]> = this.outputs$.asObservable();

  resetOutputs(): void { }

  executeBlocks(blocks: Block[]): void {
    const outputs = [];
    const output: Output = {
      blockId: 'loaddata',
      text: 'Some text here'
    };
    outputs.push(output);
    this.outputs$.next(outputs);
  }
}
