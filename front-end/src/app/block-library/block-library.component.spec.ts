import { By } from '@angular/platform-browser';
import { ComponentFixture, TestBed } from '@angular/core/testing';
import { HttpClientTestingModule } from '@angular/common/http/testing';
import { MatCardModule } from '@angular/material/card';
import { MatSnackBarModule } from '@angular/material/snack-bar';

import { BlockId } from '../block.interface';
import { BlockLibraryComponent } from './block-library.component';
import { BlockService } from '../block.service';
import { MockBlockService } from '../mock-block.service';
import { MockOutputService } from '../mock-output.service';
import { OutputService } from '../output.service';

describe('BlockLibraryComponent', () => {
  let component: BlockLibraryComponent;
  let fixture: ComponentFixture<BlockLibraryComponent>;

  beforeEach(() => {
    TestBed.configureTestingModule({
      declarations: [BlockLibraryComponent],
      imports: [
        HttpClientTestingModule,
        MatCardModule,
        MatSnackBarModule,
      ],
      providers: [
        { provide: BlockService, useClass: MockBlockService },
        { provide: OutputService, useClass: MockOutputService },
      ],
    });
    fixture = TestBed.createComponent(BlockLibraryComponent);
    component = fixture.componentInstance;
    fixture.detectChanges();
  });

  it('should create', () => {
    expect(component).toBeTruthy();
  });

  describe('Adding Blocks From Library', () => {
    it ('should call blockService.addBlock with corresponding block ID when add button is clicked', () => {
      const blockService: BlockService = TestBed.inject(BlockService);
      spyOn(blockService, 'addBlock');
      const blockIDs: BlockId[] = ['loaddata', 'basicfiltering', 'qcplots', 'qcfiltering', 'variablegenes', 'pca', 'runumap'];
      blockIDs.forEach(blockID => {
        let button = fixture.debugElement.query(By.css('#'.concat(blockID)));
        button.triggerEventHandler('click', {});
        fixture.detectChanges();
        expect(blockService.addBlock).toHaveBeenCalledWith(blockID);
      })
      expect(blockService.addBlock).toHaveBeenCalledTimes(blockIDs.length);
    });

    it ('should be available when blocks are not being executed', () => {
      component.executingBlocks = false;
      fixture.detectChanges();
      const blockIDs: string[] = ['loaddata', 'basicfiltering', 'qcplots', 'qcfiltering', 'variablegenes', 'pca', 'runumap'];
      blockIDs.forEach(blockID => {
        expect(fixture.debugElement.query(By.css('#'.concat(blockID))).nativeElement.disabled).toEqual(false);
      })
    });

    it ('should become disabled while blocks are being executed', () => {
      component.executingBlocks = false;
      fixture.detectChanges();
      const blockIDs: string[] = ['loaddata', 'basicfiltering', 'qcplots', 'qcfiltering', 'variablegenes', 'pca', 'runumap'];
      blockIDs.forEach(blockID => {
        expect(fixture.debugElement.query(By.css('#'.concat(blockID))).nativeElement.disabled).toEqual(false);
      })
      component.executingBlocks = true;
      fixture.detectChanges(); 
      blockIDs.forEach(blockID => {
        expect(fixture.debugElement.query(By.css('#'.concat(blockID))).nativeElement.disabled).toEqual(true);
      })
    });

    it ('should become available once blocks have stopped being executed', () => {
      component.executingBlocks = true;
      fixture.detectChanges();
      const blockIDs: string[] = ['loaddata', 'basicfiltering', 'qcplots', 'qcfiltering', 'variablegenes', 'pca', 'runumap'];
      blockIDs.forEach(blockID => {
        expect(fixture.debugElement.query(By.css('#'.concat(blockID))).nativeElement.disabled).toEqual(true);
      })
      component.executingBlocks = false;
      fixture.detectChanges();
      blockIDs.forEach(blockID => {
        expect(fixture.debugElement.query(By.css('#'.concat(blockID))).nativeElement.disabled).toEqual(false);
      })
    });
  });
});
